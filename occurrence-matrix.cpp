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

typedef std::map<Kmer,size_t> dictionary; // <k-mer && reverse-complement, #kmers>
typedef std::vector<Kmer> Kmers;
//typedef std::pair<size_t,size_t> cell;

// Function to create the dictionary
// assumption: kmervect has unique entries
void dictionaryCreation(dictionary &kmerdict, Kmers &kmervect)
{
    for(size_t i = 0; i<kmervect.size(); i++)
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
    std::vector<tuple<size_t,size_t,size_t>> occurrences;
    std::vector<tuple<size_t,size_t,size_t>> transtuples;

    cout << "Input k-mers file: " << argv[1] <<endl;
    cout << "Psbsim depth: " << DEPTH << endl;
    cout << "k-mer length: " << KMER_LENGTH <<endl;
    cout << "Reference genome: " << argv[2] <<endl;

    // cout << "EDLIB_EDOP_INSERT = " << EDLIB_EDOP_INSERT << endl;
    // cout << "EDLIB_EDOP_DELETE = " << EDLIB_EDOP_DELETE << endl;
    // cout << "EDLIB_EDOP_MISMATCH = " << EDLIB_EDOP_MISMATCH << endl;  

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
    
    size_t read_id = 0; // read_id needs to be global (not just per file)
    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);
        
        size_t fillstatus = 1;
        while(fillstatus) { 
            //double time1 = omp_get_wtime();
            fillstatus = pfq->fill_block(seqs, quals, upperlimit);
            size_t nreads = seqs.size();

            //double time2 = omp_get_wtime();
            //cout << "Filled " << nreads << " reads in " << time2-time1 << " seconds "<< endl; 
            
            for(size_t i=0; i<nreads; i++) 
            {
                // remember that the last valid position is length()-1
                size_t len = seqs[i].length();
                reads.push_back(seqs[i]);
                
                // skip this sequence if the length is too short
                if(len < KMER_LENGTH)
                    continue;

                for(size_t j=0; j<=len-KMER_LENGTH; j++)  
                {
                    std::string kmerstrfromfastq = seqs[i].substr(j, KMER_LENGTH);
                    Kmer mykmer(kmerstrfromfastq.c_str());
                    // remember to use only ::rep() when building kmerdict as well
                    Kmer lexsmall = mykmer.rep();      

                    auto found = kmerdict.find(lexsmall);
                    if(found != kmerdict.end()) {
                        occurrences.push_back(std::make_tuple(read_id, found->second, found->second)); // vector<tuple<read_id,kmer_id, kmer_id>
                        transtuples.push_back(std::make_tuple(found->second, read_id, found->second));
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
    //std::vector<string>().swap(seqs);   // free memory of seqs 
    std::vector<string>().swap(quals); // free memory of quals

    cout << "Total number of reads: "<< read_id << endl;
  
    CSC<size_t, size_t> spmat(occurrences, read_id, kmervect.size(), 
                            [] (size_t & c1, size_t & c2) 
                            {   if(c1 != c2) cout << "error in MergeDuplicates()" << endl;
                                return c1;
                            });
    std::cout << "spmat created with " << spmat.nnz << " nonzeros" << endl;
    std::vector<tuple<size_t,size_t,size_t>>().swap(occurrences);    // remove memory of occurences

    CSC<size_t, size_t> transpmat(transtuples, kmervect.size(), read_id, 
                            [] (size_t & c1, size_t & c2) 
                            {   if(c1 != c2) cout << "error in MergeDuplicates()" << endl;
                                return c1;
                            });
    std::cout << "transpose(spmat) created" << endl;
    std::vector<tuple<size_t,size_t,size_t>>().swap(transtuples); // remove memory of transtuples

    spmat.Sorted();
    transpmat.Sorted();

    double start = omp_get_wtime();
    CSC<size_t, size_t> tempspmat;

    HeapSpGEMM(spmat, transpmat, 
            [] (size_t & c1, size_t & c2)
            {  if(c1 != c2) 
                cout << "error in multop()" << endl;
                return 1;
            }, [] (size_t & m1, size_t & m2)
            { 
               return m1+m2;
            }, tempspmat);
    
    //cout << tempspmat.nnz << endl;
    cout << "Multiply time: " << (omp_get_wtime()-start)/60 << " min" << endl;
    //cout << "output nnz = " << tempspmat.nnz << endl;
    double start2 = omp_get_wtime();
    LocalAlignmentTest(tempspmat, reads);  // sparse mat, seq vector
    cout << "Local alignment time: " << (omp_get_wtime()-start2)/60 << " min" << endl;
    cout << "Total time: " << (omp_get_wtime()-all)/60 << " min" << endl;

    return 0;
} 
