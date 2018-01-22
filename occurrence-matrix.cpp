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
#include "mtspgemm2017/global.h"
#include "mtspgemm2017/IO.h"
#include "mtspgemm2017/multiply.h"

#define LSIZE 16000
#define KMER_LENGTH 17
#define ITERS 10
#define DEPTH 30
//#define _ALLKMER

using namespace std;

#ifdef _ALLKMER
struct spmatType_ {

    int count = 0;   /* number of shared k-mers */
    std::vector<std::pair<int,int>> vpos; /* wanna keep all the positions */
};
#else
struct spmatType_ {

    int count = 0;   /* number of shared k-mers */
    int pos[4] = {0};  /* pos1i, pos1j, pos2i, pos2j */
};
#endif

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

typedef shared_ptr<spmatType_> spmatPtr_; // pointer to spmatType_ datastruct
typedef std::map<Kmer, int> dictionary; // <k-mer && reverse-complement, #kmers>
typedef std::vector<Kmer> Kmers;

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
    if(argc < 3)
    {
        cout << "Not enough parameters. Usage: "<< endl;
        cout << "./parse kmers-list listoffastqfiles.txt" << endl;
    }

    std::ifstream filein(argv[1]);
    FILE *fastafile;
    int elem;
    char *buffer;
    string kmerstr;
    string line;
    Kmer::set_k(KMER_LENGTH);
    dictionary kmerdict;
    vector<filedata> allfiles = GetFiles(argv[2]);
    size_t upperlimit = 10000000; // in bytes
    Kmer kmerfromstr;
    Kmers kmervect;
    vector<string> seqs;
    vector<string> quals;
    vector<string> nametags;
    readVector_ reads;
    Kmers kmersfromreads;

    vector<tuple<int,int,int>> occurrences;
    vector<tuple<int,int,int>> transtuples;
    
    cout << "Input k-mer file: " << argv[1] <<endl;
    cout << "Depth: " << DEPTH << "X" << endl;
    cout << "k-mer length: " << KMER_LENGTH <<endl;
    // cout << "Sample name: PBcR-PB-ec" << endl;

    double all = omp_get_wtime();
    double kdict = omp_get_wtime();
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
    } else std::cout << "Unable to open the input file\n";
    filein.close();

    dictionaryCreation(kmerdict, kmervect);
    cout << "Reliable k-mer: " << kmerdict.size() << endl;
    cout << "k-mer dictionary creation took: " << omp_get_wtime()-kdict << "s" << endl;
    
    double parsefastq = omp_get_wtime();
    int read_id = 0; // read_id needs to be global (not just per file)
    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);
        readType_ temp;
        
        size_t fillstatus = 1;
        while(fillstatus) { 
            double time1 = omp_get_wtime();
            fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
            int nreads = seqs.size();

            // cout << "Filled " << nreads << " reads in " << time2-time1 << " seconds "<< endl; 
            
            for(int i=0; i<nreads; i++) 
            {
                // remember that the last valid position is length()-1
                int len = seqs[i].length();
                temp.nametag = nametags[i];
                temp.seq = seqs[i];
                // save reads for seeded alignment
                reads.push_back(temp);
                
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
                        occurrences.push_back(std::make_tuple(read_id, found->second, j)); // vector<tuple<read_id,kmer_id,kmerpos>>
                        transtuples.push_back(std::make_tuple(found->second, read_id, j)); // transtuples.push_back(col_id, row_id, kmerpos)
                    }
                }
                read_id++;
            }
        // cout << "processed reads in " << omp_get_wtime()-time2 << " seconds "<< endl; 
        // cout << "total number of reads processed so far is " << read_id << endl;
        } 
    delete pfq;
    }

    // don't free this vector I need this information to align sequences 
    // std::vector<string>().swap(seqs);   // free memory of seqs 
    
    std::vector<string>().swap(quals);     // free memory of quals
    cout << "Total number of reads: "<< read_id << endl;
    cout << "Parsing fastq took: " << omp_get_wtime()-parsefastq << "s" << endl;
    double matcreat = omp_get_wtime();

    #ifdef _ALLKMER
    CSC<int, int> spmat(occurrences, read_id, kmervect.size(), 
                            [] (int & p1, int & p2) 
                            {   // assume no errors in MergeDuplicates
                                // keep just the first position of that k-mer in that read
                                return p1;
                            });
    //std::cout << "spmat created with " << spmat.nnz << " nonzeros" << endl;
    std::vector<tuple<int,int,int>>().swap(occurrences);    // remove memory of occurences

    CSC<int, int> transpmat(transtuples, kmervect.size(), read_id, 
                            [] (int & p1, int & p2) 
                            {  return p1;
                            });
    //std::cout << "transpose(spmat) created" << endl;
    std::vector<tuple<int,int,int>>().swap(transtuples); // remove memory of transtuples
    #else
    CSC<int, int> spmat(occurrences, read_id, kmervect.size(), 
                            [] (int & p1, int & p2) 
                            {   // assume no errors in MergeDuplicates
                                // keep just the first position of that k-mer in that read
                                return p1;
                            });
    //std::cout << "spmat created with " << spmat.nnz << " nonzeros" << endl;
    std::vector<tuple<int,int,int>>().swap(occurrences);    // remove memory of occurences

    CSC<int, int> transpmat(transtuples, kmervect.size(), read_id, 
                            [] (int & p1, int & p2) 
                            {  return p1;
                            });
    //std::cout << "transpose(spmat) created" << endl;
    std::vector<tuple<int,int,int>>().swap(transtuples); // remove memory of transtuples
    #endif

    spmat.Sorted();
    transpmat.Sorted();
    cout << "Spmat and Spmat^T creation took: " << omp_get_wtime()-matcreat << "s" << endl;
    
    #ifdef _ALLKMER
    spmatPtr_ getvaluetype(make_shared<spmatType_>());
    HeapSpGEMM(spmat, transpmat, 
            [] (int & pi, int & pj) // n-th k-mer positions on read i and on read j 
            {   spmatPtr_ value(make_shared<spmatType_>());
                value->count = 1;
                value->vpos.push_back(make_pair(pi,pj));
                return value;
            }, 
            [] (spmatPtr_ & m1, spmatPtr_ & m2)
            {   m2->count = m1->count+m2->count;
                // insert control on independent k-mer
                m2->vpos.insert(m2->vpos.end(), m1->vpos.begin(), m1->vpos.end());
                return m2;
            }, reads, getvaluetype);
    #else
    spmatPtr_ getvaluetype(make_shared<spmatType_>());
    HeapSpGEMM(spmat, transpmat, 
            [] (int & pi, int & pj) // n-th k-mer positions on read i and on read j 
            {   spmatPtr_ value(make_shared<spmatType_>());
                value->count = 1;
                value->pos[0] = pi; // row
                value->pos[1] = pj; // col
                return value;
            }, 
            [] (spmatPtr_ & m1, spmatPtr_ & m2)
            {   m2->count = m1->count+m2->count;
                // m1->pos[0] = m1->pos[0]; // row 
                // m1->pos[1] = m1->pos[1]; // col
                m2->pos[2] = m1->pos[0]; // row 
                m2->pos[3] = m1->pos[1]; // col
                return m2;
            }, reads, getvaluetype);
    #endif 
    
    cout << "Total running time: " << omp_get_wtime()-all << "s" << endl;

    return 0;
} 