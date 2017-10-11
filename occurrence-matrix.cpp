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
#include "mtspgemm2017/IO.h"
#include "mtspgemm2017/global.h"
#include "mtspgemm2017/multiply.h"

#define LSIZE 16000
#define KMER_LENGTH 17
#define ITERS 10
#define DEPTH 25
#define _SEEDED
//#define _MULPTR

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

#ifdef _MULPTR
typedef std::pair<int, vector<int>> cellspmat;             // pair<kmer_id_j, vector<posix_in_read_i>>
typedef std::vector<pair<int, pair<int, int>>> multcell;   // map<kmer_id, vector<posix_in_read_i, posix_in_read_j>>
typedef shared_ptr<multcell> mult_ptr;                     // pointer to multcell
#endif

#ifdef _SEEDED
struct spmatype {

    int count = 0;   /* number of shared k-mers */
    int pos[4] = {0};  /* pos1i, pos1j, pos2i, pos2j, pos3i, pos3j */
};
typedef shared_ptr<spmatype> spmat_ptr; // pointer to spmatype datastruct
struct readsinfo {

    std::string nametag;   
    std::string seq; 

};
#endif

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
    std::vector<string> quals;
    std::vector<string> nametags;
    std::vector<readsinfo> reads;
    int rangeStart;
    readsinfo temp;
    Kmers kmersfromreads;

    #ifdef _MULPTR
    std::vector<tuple<int,int,cellspmat>> occurrences;
    std::vector<tuple<int,int,cellspmat>> transtuples;
    #endif
    #ifdef _SEEDED
    std::vector<tuple<int,int,int>> occurrences;
    std::vector<tuple<int,int,int>> transtuples;
    #else
    std::vector<tuple<int,int,std::pair<int,int>>> occurrences; // I could need just int also in this light version to keep track of the k-mer position in the read
    std::vector<tuple<int,int,std::pair<int,int>>> transtuples;
    #endif
    
    cout << "input k-mers file: " << argv[1] <<endl;
    cout << "pbsim depth: " << DEPTH << endl;
    cout << "k-mer length: " << KMER_LENGTH <<endl;
    cout << "reference genome: " << argv[2] <<endl;

    #ifdef _MULPTR
    cout << "*** MULPTR version ***" << endl;
    #endif

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
    cout << "reliable k-mers = " << kmerdict.size() << endl;
    
    int read_id = 0; // read_id needs to be global (not just per file)
    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);
        
        size_t fillstatus = 1;
        while(fillstatus) { 
            //double time1 = omp_get_wtime();
            fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
            int nreads = seqs.size();

            //double time2 = omp_get_wtime();
            //cout << "Filled " << nreads << " reads in " << time2-time1 << " seconds "<< endl; 
            
            for(int i=0; i<nreads; i++) 
            {
                // remember that the last valid position is length()-1
                int len = seqs[i].length();
                temp.nametag = nametags[i];
                temp.seq = seqs[i];
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
                        #ifdef _MULPTR
                        occurrences.push_back(std::make_tuple(read_id, found->second, make_pair(found->second, vector<int>(1,j)))); // vector<tuple<read_id,kmer_id,pair<kmer_id,pos_in_read>>
                        transtuples.push_back(std::make_tuple(found->second, read_id, make_pair(found->second, vector<int>(1,j)))); // transtuples.push_back(col_id, row_id, value)
                        #endif
                        #ifdef _SEEDED
                        occurrences.push_back(std::make_tuple(read_id, found->second, j)); // vector<tuple<read_id,kmer_id,kmerpos>>
                        transtuples.push_back(std::make_tuple(found->second, read_id, j)); // transtuples.push_back(col_id, row_id, kmerpos)
                        #else
                        occurrences.push_back(std::make_tuple(read_id, found->second, make_pair(found->second, j))); // vector<tuple<read_id,kmer_id,kmer_id + kmerpos>
                        transtuples.push_back(std::make_tuple(found->second, read_id, make_pair(found->second, j)));
                        #endif
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
    cout << "total number of reads: "<< read_id << endl;

    #ifdef _MULPTR
    CSC<int, cellspmat> spmat(occurrences, read_id, kmervect.size(), 
                            [] (cellspmat & c1, cellspmat & c2) 
                            {   if(c1.first != c2.first) cout << "error in MergeDuplicates()" << endl;
                                vector<int> merged;
                                merge(c1.second.begin(), c1.second.end(), c2.second.begin(), c2.second.end(), back_inserter(merged));
                                return make_pair(c1.first, merged);
                            });
    std::vector<tuple<int,int,cellspmat>>().swap(occurrences);    // remove memory of occurences

    CSC<int, cellspmat> transpmat(transtuples, kmervect.size(), read_id, 
                            [] (cellspmat & c1, cellspmat & c2) 
                            {   if(c1.first != c2.first) cout << "error in MergeDuplicates()" << endl;
                                vector<int> merged;
                                merge(c1.second.begin(), c1.second.end(), c2.second.begin(), c2.second.end(), back_inserter(merged));
                                return make_pair(c1.first, merged);
                            });
    std::vector<tuple<int,int,cellspmat>>().swap(transtuples); // remove memory of transtuples
    #endif
    #ifdef _SEEDED
    CSC<int, int> spmat(occurrences, read_id, kmervect.size(), 
                            [] (int & p1, int & p2) 
                            {   // I assume that there's no errors in MergeDuplicates (to fix)
                                // I keep just the first position of that k-mer in that read
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
    CSC<int, std::pair<int, int>> spmat(occurrences, read_id, kmervect.size(), 
                            [] (std::pair<int,int> & c1, std::pair<int,int> & c2) 
                            {   if(c1.first != c2.first) cout << "error in MergeDuplicates()" << endl;
                                std::pair<int, int> value;
                                value = std::make_pair(c1.first, c1.second);
                                return value;
                            });
    //std::cout << "spmat created with " << spmat.nnz << " nonzeros" << endl;
    std::vector<tuple<int,int,std::pair<int, int>>>().swap(occurrences);    // remove memory of occurences

    CSC<int, std::pair<int, int>> transpmat(transtuples, kmervect.size(), read_id, 
                            [] (std::pair<int,int> & c1, std::pair<int,int> & c2) 
                            {   if(c1.first != c2.first) cout << "error in MergeDuplicates()" << endl;
                                std::pair<int, int> value;
                                value = std::make_pair(c1.first, c1.second);
                                return value;
                            });
    //std::cout << "transpose(spmat) created" << endl;
    std::vector<tuple<int,int,std::pair<int, int>>>().swap(transtuples); // remove memory of transtuples
    #endif

    spmat.Sorted();
    transpmat.Sorted();

    double start = omp_get_wtime();
    
    #ifdef _MULPTR
    mult_ptr getvaluetype(make_shared<multcell>());
    HeapSpGEMM(spmat, transpmat, 
        [] (cellspmat & c1, cellspmat & c2)
        {  if(c1.first != c2.first) cout << "error in multop()" << endl; 
                mult_ptr value(make_shared<multcell>()); // only one allocation
                for(int i=0; i<c1.second.size(); ++i) {
                    for(int j=0; j<c2.second.size(); ++j) {
                        pair<int, int> temp = make_pair(c1.second[i], c2.second[j]);
                        value->push_back(make_pair(c1.first, temp));
                    }
                }
                return value;
        }, 
        [] (mult_ptr & h, mult_ptr & m)
            {   m->insert(m->end(), h->begin(), h->end());
            return m;
            }, reads, getvaluetype);
    #endif
    #ifdef _SEEDED
    spmat_ptr getvaluetype(make_shared<spmatype>());
    HeapSpGEMM(spmat, transpmat, 
            [] (int & pi, int & pj) // n-th k-mer positions on read i and on read j 
            {   spmat_ptr value(make_shared<spmatype>());
                value->count = 1;
                value->pos[0] = pi; // row
                value->pos[1] = pj; // col
                return value;
            }, 
            [] (spmat_ptr & m1, spmat_ptr & m2)
            {   spmat_ptr value(make_shared<spmatype>());
                value->count = m1->count+m2->count;
                value->pos[0] = m1->pos[0]; // row 
                value->pos[1] = m1->pos[1]; // col
                value->pos[2] = m2->pos[0]; // row 
                value->pos[3] = m2->pos[1]; // col
                return value;
            }, reads, getvaluetype);
    #else 
    std::pair<int, std::pair<int, int>> getvaluetype;
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
            }, reads, getvaluetype);
    #endif
    
    cout << "multiply (local alignment included) time: " << omp_get_wtime()-start << " sec" << endl;
    // std::ifstream filename("spmat.csv");
    // getMetrics(filename);  
    cout << "total time: " << omp_get_wtime()-all << " sec" << endl;

    return 0;
} 
