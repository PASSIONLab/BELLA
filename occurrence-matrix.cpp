#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <array>
#include <tuple>
#include <queue>
#include <memory>
#include <stack>
#include <functional>
#include <cstring>
#include <string.h>
#include <cassert>
#include <ios>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <map>
#include <omp.h>
#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include "kmercode/Buffer.h"
#include "kmercode/common.h"
#include "kmercode/fq_reader.h"
#include "kmercode/ParallelFASTQ.h"

#include "mtspgemm2017/utility.h"
#include "mtspgemm2017/CSC.h"
#include "mtspgemm2017/CSR.h"
#include "edlib/edlib/include/edlib.h"
#include "mtspgemm2017/IO.h"
#include "mtspgemm2017/global.h"
#include "mtspgemm2017/metric.h"
#include "mtspgemm2017/multiply.h"

#define LSIZE 16000
#define KMER_LENGTH 17
#define ITERS 10
#define DEPTH 30

using namespace std;

struct filedata {

    char filename[MAX_FILE_PATH];
    size_t filesize;
};

std::vector<filedata>  GetFiles(char *filename) {
    int64_t totalsize = 0;
    int numfiles = 0;
    std::vector<filedata> filesview;
    
    filedata fdata;
    ifstream allfiles(filename);
    if(!allfiles.is_open()) {
        cerr << "could not open " << filename << endl;
        exit(1);
    }
    allfiles.getline(fdata.filename,MAX_FILE_PATH);
    while(!allfiles.eof())
    {
        struct stat st;
        stat(fdata.filename, &st);
        fdata.filesize = st.st_size;
        
        filesview.push_back(fdata);
        cout << filesview.back().filename << " : " << filesview.back().filesize / (1024*1024) << " MB" << endl;
        allfiles.getline(fdata.filename,MAX_FILE_PATH);
        totalsize += fdata.filesize;
        numfiles++;
    }
    return filesview;
}

typedef std::map<Kmer, int> dictionary; // <k-mer && reverse-complement, #kmers>
typedef std::vector<Kmer> Kmers;
//typedef std::pair<int, int> cell; // pair<kmer_id_j, vector<posix_in_read_i>>
//typedef std::vector<pair<size_t, pair<size_t, size_t>>> multcell; // map<kmer_id, vector<posix_in_read_i, posix_in_read_j>>
//typedef shared_ptr<multcell> mult_ptr; // pointer to multcell

// Function to create the dictionary
// assumption: kmervect has unique entries
void dictionaryCreation(dictionary &kmerdict, Kmers &kmervect)
{
    for(int i = 0; i<kmervect.size(); i++)
    {
        kmerdict.insert(make_pair(kmervect[i].rep(), i));
    }
    
}

int main (int argc, char* argv[]) {
    if(argc < 4)
    {
        cout << "not enough parameters... usage: "<< endl;
        cout << "./parse kmers-list reference-genome.fna listoffastqfiles.txt" << endl;
    }

    std::ifstream filein(argv[1]);
    FILE *fastafile;
    int elem;
    char *buffer;
    std::string kmerstr;
    std::string line;
    Kmer::set_k(KMER_LENGTH);
    dictionary kmerdict;
    std::vector<filedata> allfiles = GetFiles(argv[3]);
    size_t upperlimit = 10000000; // in bytes
    Kmer kmerfromstr;
    Kmers kmervect;
    std::vector<string> seqs;
    std::vector<string> reads;
    std::vector<string> quals;
    int rangeStart;
    Kmers kmersfromreads;
    std::vector<tuple<int,int,std::pair<int, int>>> occurrences;
    std::vector<tuple<int,int,std::pair<int, int>>> transtuples;
    
    cout << "Input k-mers file: " << argv[1] <<endl;
    cout << "Psbsim depth: " << DEPTH << endl;
    cout << "k-mer length: " << KMER_LENGTH <<endl;
    cout << "Reference genome: " << argv[2] <<endl;

    double all = omp_get_wtime();
    if(filein.is_open()) {
            while(getline(filein, line)) {
                if(line.length() == 0)
                    break;

                string substring = line.substr(1);
                elem = stoi(substring);
                getline(filein, kmerstr);   
                kmerfromstr.set_kmer(kmerstr.c_str());
                kmervect.push_back(kmerfromstr);                                
            }
    } else std::cout << "unable to open the input file\n";
    filein.close();

    dictionaryCreation(kmerdict, kmervect);
    cout << "Reliable k-mers = " << kmerdict.size() << endl;
    
    int read_id = 0; // read_id needs to be global (not just per file)
    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);
        
        size_t fillstatus = 1;
        while(fillstatus) { 
            //double time1 = omp_get_wtime();
            fillstatus = pfq->fill_block(seqs, quals, upperlimit);
            int nreads = seqs.size();

            //double time2 = omp_get_wtime();
            //cout << "Filled " << nreads << " reads in " << time2-time1 << " seconds "<< endl; 
            
            for(int i=0; i<nreads; i++) 
            {
                // remember that the last valid position is length()-1
                int len = seqs[i].length();
                reads.push_back(seqs[i]);
                
                // skip this sequence if the length is too short
                if(len < KMER_LENGTH)
                    continue;

                for(int j=0; j<=len-KMER_LENGTH; j++)  
                {
                    std::string kmerstrfromfastq = seqs[i].substr(j, KMER_LENGTH);
                    Kmer mykmer(kmerstrfromfastq.c_str());
                    // remember to use only ::rep() when building kmerdict as well
                    Kmer lexsmall = mykmer.rep();      

                    auto found = kmerdict.find(lexsmall);
                    if(found != kmerdict.end()) {
                        occurrences.push_back(std::make_tuple(read_id, found->second, make_pair(found->second, j))); // vector<tuple<read_id,kmer_id, kmer_id>
                        transtuples.push_back(std::make_tuple(found->second, read_id, make_pair(found->second, j)));
                    }
                }
                read_id++;
            }
        //cout << "processed reads in " << omp_get_wtime()-time2 << " seconds "<< endl; 
        //cout << "total number of reads processed so far is " << read_id << endl;
        } 
    delete pfq;
    }
    // don't free this vector I need this information to align sequences 
    // std::vector<string>().swap(seqs);   // free memory of seqs 
    std::vector<string>().swap(quals);     // free memory of quals

    cout << "Total number of reads: "<< read_id << endl;
  
    CSC<int, std::pair<int, int>> spmat(occurrences, read_id, kmervect.size(), 
                            [] (std::pair<int,int> & c1, std::pair<int,int> & c2) 
                            {   if(c1.first != c2.first) cout << "error in MergeDuplicates()" << endl;
                                std::pair<int, int> value;
                                value = std::make_pair(c1.first, c1.second);
                                return value;
                            });
    std::cout << "spmat created with " << spmat.nnz << " nonzeros" << endl;
    std::vector<tuple<int,int,std::pair<int, int>>>().swap(occurrences);    // remove memory of occurences

    CSC<int, std::pair<int, int>> transpmat(transtuples, kmervect.size(), read_id, 
                            [] (std::pair<int,int> & c1, std::pair<int,int> & c2) 
                            {   if(c1.first != c2.first) cout << "error in MergeDuplicates()" << endl;
                                std::pair<int, int> value;
                                value = std::make_pair(c1.first, c1.second);
                                return value;
                            });
    std::cout << "transpose(spmat) created" << endl;
    std::vector<tuple<int,int,std::pair<int, int>>>().swap(transtuples); // remove memory of transtuples

    spmat.Sorted();
    transpmat.Sorted();

    double start = omp_get_wtime();
    
    /*CSC<size_t, mult_ptr> tempspmat;
    cout << "before multiply" <<endl;

    HeapSpGEMM(spmat, transpmat, 
        [] (cellspmat & c1, cellspmat & c2)
        {  if(c1.first != c2.first) cout << "error in multop()" << endl; 
                mult_ptr value(make_shared<multcell>()); // only one allocation
                for(int i=0; i<c1.second.size(); ++i) {
                    for(int j=0; j<c2.second.size(); ++j) {
                        pair<size_t, size_t> temp = make_pair(c1.second[i], c2.second[j]);
                        value->push_back(make_pair(c1.first, temp));
                    }
                }
                return value;
        }, 
        [] (mult_ptr & h, mult_ptr & m)
            {   m->insert(m->end(), h->begin(), h->end());
            return m;
            }, tempspmat);
    
    cout << "multiply took " << omp_get_wtime()-start << " seconds" << endl;
    double start2 = omp_get_wtime();
    DetectOverlap(tempspmat); // function to remove reads pairs that don't show evidence of potential overlap and to compute the recall 
    cout << "Filter time " << omp_get_wtime()-start2 << " seconds" << endl;
    cout << "Total time " << omp_get_wtime()-all << " seconds" << endl; */

    // CSC<size_t, size_t> tempspmat; 
    HeapSpGEMM(spmat, transpmat, 
            [] (std::pair<int,int> & c1, std::pair<int,int> & c2)
            {   if(c1.first != c2.first) cout << "error in multop()" << endl;
                std::pair<int, int> pos;
                pos = std::make_pair(c1.second, c2.second);
                std::pair<int, std::pair<int, int>> value;
                value = std::make_pair(1, pos);
                return value;
            }, 
            [] (std::pair<int, std::pair<int, int>> & m1, std::pair<int, std::pair<int, int>> & m2)
            {   std::pair<int, std::pair<int, int>> value;
                value = std::make_pair(m1.first+m2.first, std::make_pair(m2.second.first, m2.second.second));
                return value;
            }, reads);
    
    cout << "Multiply time: " << (omp_get_wtime()-start)/60 << " min" << endl;
    std::ifstream filename("spmat.csv");
    GetMetrics(filename);  
    cout << "Total time: " << (omp_get_wtime()-all)/60 << " min" << endl;
    return 0;
} 
