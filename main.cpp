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
#include <math.h>
#include <cassert>
#include <ios>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <map>
#include <unordered_map>
#include <omp.h>

#include "libcuckoo/cuckoohash_map.hh"
#include "kmercount.h"

#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include "kmercode/Buffer.h"
#include "kmercode/common.h"
#include "kmercode/fq_reader.h"
#include "kmercode/ParallelFASTQ.h"
#include "kmercode/bound.hpp"

#include "mtspgemm2017/utility.h"
#include "mtspgemm2017/CSC.h"
#include "mtspgemm2017/CSR.h"
#include "mtspgemm2017/global.h"
#include "mtspgemm2017/IO.h"
#include "mtspgemm2017/overlapping.h"
#include "mtspgemm2017/align.h"

#define LSIZE 16000
#define ITERS 10
#define PRINT
//#define JELLYFISH
//#define ADAPTIVE TODO: ifdef adaptive version otherwise default

using namespace std;

int main (int argc, char *argv[]) {

    //
    // Program name and purpose
    //
    cout << "\nBELLA - Long Read Aligner for De Novo Genome Assembly\n" << endl;
    //
    // Setup the input files
    //
    option_t *optList, *thisOpt;
    // Get list of command line options and their arguments 
    optList = NULL;
    optList = GetOptList(argc, argv, (char*)"f:i:o:d:hk:a:ze:p:w:m:y:g:s:r:");
   
    bool skip_algnmnt_krnl = false;     // Do not align (z)
    char *kmer_file = NULL;             // Reliable k-mer file from Jellyfish
    char *all_inputs_fofn = NULL;       // List of fastqs (i)
    char *out_file = NULL;              // output filename (o)
    int kmer_len = 17;                  // default k-mer length (k)
    int algnmnt_thr = 50;               // default alignment score threshold (a)
    int xdrop = 7;                      // default alignment x-drop factor (p)
    double erate = 0.15;                // default error rate (e)
    int depth = 0;                      // depth/coverage required (d)
    int epsilon = 300;                  // epsilon parameter for alignment on edges TODO: explain (w)
    vector<int> scores = {1,-1,-1,-7,9};  // mat (m), mis (s), gex (y), gop (g), probability 1:9 substitution:indel ratio (r)


    if(optList == NULL)
    {
        cout << "BELLA execution terminated: not enough parameters or invalid option" << endl;
        cout << "Run with -h to print out the command line options\n" << endl;
        return 0;
    }

    while (optList!=NULL) {
        thisOpt = optList;
        optList = optList->next;
        switch (thisOpt->option) {
            case 'f': {
#ifdef JELLYFISH
                if(thisOpt->argument == NULL)
                {
                    cout << "BELLA execution terminated: -f requires an argument" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
#endif
                kmer_file = strdup(thisOpt->argument);
                break;
            }
            case 'i': {
                if(thisOpt->argument == NULL)
                {
                    cout << "BELLA execution terminated: -i requires an argument" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                all_inputs_fofn = strdup(thisOpt->argument);
                break;
            }
            case 'o': {
                if(thisOpt->argument == NULL)
                {
                    cout << "BELLA execution terminated: -o requires an argument" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                char* line1 = strdup(thisOpt->argument);
                char* line2 = strdup(".out");
                size_t len1 = strlen(line1);
                size_t len2 = strlen(line2);

                out_file = (char*)malloc(len1 + len2 + 1);
                if (!out_file) abort();

                memcpy(out_file, line1, len1);
                memcpy(out_file + len1, line2, len2);
                out_file[len1 + len2] = '\0';

                delete line1;
                delete line2;

                break;
            }
            case 'd': {
                if(thisOpt->argument == NULL)
                {
                    cout << "BELLA execution terminated: -d requires an argument" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                depth = atoi(thisOpt->argument);  
                break;
            }
            case 'z': skip_algnmnt_krnl = true; break; // TO DO: add skip alignment
            case 'k': {
                kmer_len = atoi(thisOpt->argument);
                break;
            }
            case 'e': {
                erate = strtod(thisOpt->argument, NULL);
                break;
            }
            case 'a': {
                algnmnt_thr = atoi(thisOpt->argument);
                break;
            }
            case 'p': {
                xdrop = atoi(thisOpt->argument);
                break;
            }
            case 'w': {
                epsilon = atoi(thisOpt->argument);
                break;
            }
            case 'm': {
                if(thisOpt->argument < 0)
                {
                    cout << "BELLA execution terminated: -m requires a positive integer" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                scores[0] = atoi(thisOpt->argument);
                break;
            }
            case 's': {
                if(thisOpt->argument > 0)
                {
                    cout << "BELLA execution terminated: -s requires a negative integer" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                scores[1] = atoi(thisOpt->argument);
                break;
            }
            case 'y': {
                if(thisOpt->argument > 0)
                {
                    cout << "BELLA execution terminated: -y requires a negative integer" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                scores[2] = atoi(thisOpt->argument);
                break;
            }
            case 'g': {
                if(thisOpt->argument > 0)
                {
                    cout << "BELLA execution terminated: -g requires a negative integer" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                scores[4] = atoi(thisOpt->argument);
                break;
            }
            case 'r': {
                if(thisOpt->argument < 0)
                {
                    cout << "BELLA execution terminated: -r requires a positive integer" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                scores[5] = atoi(thisOpt->argument);
                break;
            }
            case 'h': {
                cout << "Usage:\n" << endl;
                cout << " -f : k-mer list from Jellyfish (required if Jellyfish k-mer counting is used)" << endl; // the reliable k-mers are selected by bella
                cout << " -i : list of fastq(s) (required)" << endl;
                cout << " -o : output filename (required)" << endl;
                cout << " -d : depth/coverage (required)" << endl; // TO DO: add depth estimation
                cout << " -k : k-mer length [17]" << endl;
                cout << " -a : alignment score threshold [50]" << endl;
                cout << " -p : alignment x-drop factor [3]" << endl;
                cout << " -e : error rate [0.15]" << endl;
                cout << " -z : skip the alignment [false]\n" << endl;
                cout << " -w : epsilon parameter for alignment on edges [300]\n" << endl;
                cout << " -m : match penalty scoring matrix [1]\n" << endl;
                cout << " -s : mismatch penalty scoring matrix [-1]\n" << endl;
                cout << " -y : gap extension penalty scoring matrix [-1]\n" << endl;
                cout << " -g : gap opening penalty scoring matrix [-7]\n" << endl;
                cout << " -r : substitution:indel probability ratio [1:9=sub:gap]\n" << endl;

                FreeOptList(thisOpt); // Done with this list, free it
                return 0;
            }
        }
    }

#ifdef JELLYFISH
    if(kmer_file == NULL || all_inputs_fofn == NULL || out_file == NULL || depth == 0)
    {
        cout << "BELLA execution terminated: missing arguments" << endl;
        cout << "Run with -h to print out the command line options\n" << endl;
        return 0;
    }
#else
    if(all_inputs_fofn == NULL || out_file == NULL || depth == 0)
    {
        cout << "BELLA execution terminated: missing arguments" << endl;
        cout << "Run with -h to print out the command line options\n" << endl;
        return 0;
    }
#endif

    free(optList);
    free(thisOpt);
    //
    // Declarations 
    //
    vector<filedata> allfiles = GetFiles(all_inputs_fofn);
    FILE *fastafile;
    int lower; // reliable range lower bound (fixed)
    int upper;     // reliable range upper bound
    char *buffer;
    Kmer::set_k(kmer_len);
    size_t upperlimit = 10000000; // in bytes
    Kmers kmervect;
    vector<string> seqs;
    vector<string> quals;
    vector<string> nametags;
    readVector_ reads;
    Kmers kmersfromreads;
    vector<tuple<int,int,int>> occurrences;
    vector<tuple<int,int,int>> transtuples;
    // 
    // File and setting used
    //
#ifdef PRINT
#ifdef JELLYFISH
    cout << "K-mer counting: Jellyfish" << endl;
    cout << "Input k-mer file: " << kmer_file << endl;
#endif
    cout << "K-mer counting: BELLA" << endl;
    //cout << "minimum number of shared k-mer: " << n << endl;
    cout << "Epsilon: " << epsilon << endl;
    cout << "Output filename: " << out_file << endl;
    cout << "K-mer length: " << kmer_len << endl;
    if(skip_algnmnt_krnl)
        cout << "Compute alignment: false" << endl;
    else cout << "Compute alignment: true" << endl;
    cout << "X-drop: " << xdrop << endl;
    cout << "Default alignment score threshold: " << algnmnt_thr << endl;
    cout << "Depth: " << depth << "X" << endl;
    cout << "Error rate: " << erate << endl;
    cout << "Scoring matrix: (" << scores[0] << "," << scores[1] << "," << scores[2] << "," << scores[3] << ")" << endl;
    cout << "Ratio substitution:indel = 1:" << scores[4] << endl;
#endif
    //
    // Reliable bounds computation
    //
    double all = omp_get_wtime();
    lower = computeLower(depth,erate,kmer_len);
    upper = computeUpper(depth,erate,kmer_len);
    double A = adaptiveSlope(erate, (float)xdrop, scores);
#ifdef PRINT
    cout << "Reliable lower bound: " << lower << endl;
    cout << "Reliable upper bound: " << upper << "\n" << endl;
    cout << "Constant of adaptive threshold: " << A << endl;
#endif
    //
    // Reliable k-mer file parsing and k-mer dictionary creation
    //
    dictionary_t countsreliable;
#ifdef JELLYFISH
    JellyFishCount(kmer_file, countsreliable, lower, upper);
#else
    cout << "\nRunning with up to " << MAXTHREADS << " threads" << endl;
    DeNovoCount(allfiles, countsreliable, lower, upper, kmer_len, upperlimit);
#endif
    //
    // Fastq(s) parsing
    //
    double parsefastq = omp_get_wtime();
    int read_id = 0; // read_id needs to be global (not just per file)
    cout << "\nRunning with up to " << MAXTHREADS << " threads" << endl;

        vector < vector<tuple<int,int,int>> > alloccurrences(MAXTHREADS);   
        vector < vector<tuple<int,int,int>> > alltranstuples(MAXTHREADS);   
        vector < readVector_ > allreads(MAXTHREADS);

        double time1 = omp_get_wtime();
        for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);

        size_t fillstatus = 1;
        while(fillstatus) { 
            fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
            int nreads = seqs.size();
       
                #pragma omp parallel for
                for(int i=0; i<nreads; i++) 
                {
                    // remember that the last valid position is length()-1
                    int len = seqs[i].length();

                    readType_ temp;
                    temp.nametag = nametags[i];
                    temp.seq = seqs[i]; // save reads for seeded alignment
                    temp.readid = read_id+i;

		    allreads[MYTHREAD].push_back(temp);
                    
                    for(int j=0; j<=len-kmer_len; j++)  
                    {
                        std::string kmerstrfromfastq = seqs[i].substr(j, kmer_len);
                        Kmer mykmer(kmerstrfromfastq.c_str());
                        // remember to use only ::rep() when building kmerdict as well
                        Kmer lexsmall = mykmer.rep();      

                        int idx; // kmer_id
			auto found = countsreliable.find(lexsmall,idx);
			if(found)
                        {
                            alloccurrences[MYTHREAD].emplace_back(std::make_tuple(read_id+i,idx,j)); // vector<tuple<read_id,kmer_id,kmerpos>>
                            alltranstuples[MYTHREAD].emplace_back(std::make_tuple(idx,read_id+i,j)); // transtuples.push_back(col_id,row_id,kmerpos)
                        }
                    }
                }   // for(int i=0; i<nreads; i++)
                    //cout << "total number of reads processed so far is " << read_id << endl;
            read_id += nreads;
        } //while(fillstatus) 
    delete pfq;
    } // for all files

    size_t readcount = 0;
    size_t tuplecount = 0;
    for(int t=0; t<MAXTHREADS; ++t)
    {
        readcount += allreads[t].size();
        tuplecount += alloccurrences[t].size();
    }
    reads.resize(readcount);
    occurrences.resize(tuplecount);
    transtuples.resize(tuplecount);

    size_t readssofar = 0;
    size_t tuplesofar = 0;
    for(int t=0; t<MAXTHREADS; ++t)
    {
        copy(allreads[t].begin(), allreads[t].end(), reads.begin()+readssofar);
        readssofar += allreads[t].size();

        copy(alloccurrences[t].begin(), alloccurrences[t].end(), occurrences.begin() + tuplesofar);
        copy(alltranstuples[t].begin(), alltranstuples[t].end(), transtuples.begin() + tuplesofar);
        tuplesofar += alloccurrences[t].size();
    }

    std::sort(reads.begin(), reads.end());   // bool operator in global.h: sort by readid
    std::vector<string>().swap(seqs);        // free memory of seqs  
    std::vector<string>().swap(quals);       // free memory of quals
    cout << "\nTotal number of reads: "<< read_id << endl;
    cout << "fastq(s) parsing fastq took: " << omp_get_wtime()-parsefastq << "s" << endl;
    double matcreat = omp_get_wtime();
    //
    // Sparse matrices construction
    //
    int nkmer = countsreliable.size();
    CSC<int, int> spmat(occurrences, read_id, nkmer, 
                            [] (int & p1, int & p2) 
                            {   // assume no errors in MergeDuplicates
                                // keep just the first position of that k-mer in that read
                                return p1;
                            });
    std::vector<tuple<int,int,int>>().swap(occurrences);    // remove memory of occurences

    CSC<int, int> transpmat(transtuples, nkmer, read_id, 
                            [] (int & p1, int & p2) 
                            {  return p1;
                            });
    std::vector<tuple<int,int,int>>().swap(transtuples); // remove memory of transtuples

    spmat.Sorted();
    transpmat.Sorted();

#ifdef PRINT
    cout << "spmat and spmat^T creation took: " << omp_get_wtime()-matcreat << "s" << endl;
#endif
    //
    // Overlap detection (sparse matrix multiplication) and seed-and-extend alignment
    //
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
                m2->pos[2] = m1->pos[0]; // row 
                m2->pos[3] = m1->pos[1]; // col
                return m2;
            }, reads, getvaluetype, kmer_len, xdrop, algnmnt_thr, out_file, skip_algnmnt_krnl, scores, A); 

    cout << "total running time: " << omp_get_wtime()-all << "s\n" << endl;
    return 0;
} 
