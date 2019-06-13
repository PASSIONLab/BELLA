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
#include "mtspgemm2017/common.h"
#include "mtspgemm2017/IO.h"
#include "mtspgemm2017/overlapping.h"
#include "mtspgemm2017/align.h"

#define LSIZE 16000
#define ITERS 10
#define PRINT
//#define JELLYFISH

double deltaChernoff = 0.1;            

using namespace std;



void testKmer() {
  // TODO: remove this, this is just to test that the Kmer class works right with different lengths, etc
  
  cout << "Testing Kmer length" << endl;
  // create kmer of the default length (17)
    // with the default constructor
  Kmer def();
    // with the string constructor
  Kmer str("AAAAAAAAAAAAAAAAA", 17);
    // with the copy constructor
  // Kmer copy(&def);
    // with the assignment op
  Kmer assign = str;
  cout << "Length k passed" << endl;
  cout << "HOPC: " << str.getHOPC() << endl;
  
  // create kmer of length 0
    // with the string constructor
  Kmer empty("", 0);
    // with the copy and assignment constructors
  // Kmer empty_copy(&empty);
  Kmer empty_assign = empty;
  cout << "Empty passed" << endl;
  cout << "HOPC: " << empty.getHOPC() << endl;
  
  // create kmer of maximum length
    // with the string constructor
    // TODO: figure out how to do this

  // create kmer of length 5
  Kmer str5("ACTTG", 5);
  // Kmer copy5(&str5);
  Kmer assign5 = str5;
  cout << "Length 5 passed" << endl;
  cout << "HOPC: " << str5.getHOPC() << endl;

  // create kmer of length 22
  Kmer str20("ACTGGGGGGGTTTTAAAACCCC", 22);
  // Kmer copy20(&str20);
  Kmer assign20 = str20;
  cout << "Length 22 passed" << endl;
  cout << "HOPC: " << str20.getHOPC() << endl;

  // perform comparison operations on these kmers


  // decide if we need to be able to change the length of existing Kmers
}

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
    // Follow an option with a colon to indicate that it requires an argument.

    optList = NULL;
    optList = GetOptList(argc, argv, (char*)"f:i:o:d:hk:Ka:ze:x:w:nc:m:r:p:u");
   

    char *kmer_file = NULL;                 // Reliable k-mer file from Jellyfish
    char *all_inputs_fofn = NULL;           // List of fastqs (i)
    char *out_file = NULL;                  // output filename (o)
    int kmer_len = 17;                      // default k-mer length (k)
    int xdrop = 7;                          // default alignment x-drop factor (x)
    double erate = 0.15;                    // default error rate (e) 
    int depth = 0;                          // depth/coverage required (d)

    BELLApars b_parameters;

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
            case 'p': b_parameters.outputPaf = true; break; // PAF format
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
            case 'z': b_parameters.skipAlignment = true; break;
            case 'n': b_parameters.alignEnd = true; break; 
            case 'k': {
                kmer_len = atoi(thisOpt->argument);
                break;
            }
            case 'r': {
                b_parameters.kmerRift = atoi(thisOpt->argument);
                break;
            }
            case 'e': {
                b_parameters.skipEstimate = true;
                erate = strtod(thisOpt->argument, NULL);
                break;
            }
            case 'a': {
                b_parameters.defaultThr = atoi(thisOpt->argument);
                b_parameters.adapThr = false;
                break;
            }
            case 'x': {
                xdrop = atoi(thisOpt->argument);
                break;
            } 
            case 'w': {
                b_parameters.relaxMargin = atoi(thisOpt->argument);
                break;
            }
            case 'K': {
                b_parameters.allKmer = true;
                break;
            }
            case 'm': {
                b_parameters.totalMemory = stod(thisOpt->argument);
                b_parameters.userDefMem = true;
                cout << "User defined memory set to " << b_parameters.totalMemory << " MB " << endl;
                break;
            }
            case 'c': {
                if(stod(thisOpt->argument) > 1.0 || stod(thisOpt->argument) < 0.0)
                {
                    cout << "BELLA execution terminated: -c requires a value in [0,1]" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                b_parameters.deltaChernoff = stod(thisOpt->argument);
                break;
            }
            case 'u': {
              b_parameters.useHOPC = true;
              break;
            }
            case 'h': {
                cout << "Usage:\n" << endl;
                cout << " -f : k-mer list from Jellyfish (required if Jellyfish k-mer counting is used)" << endl; // Reliable k-mers are selected by BELLA
                cout << " -i : list of fastq(s) (required)" << endl;
                cout << " -o : output filename (required)" << endl;
                cout << " -d : depth/coverage (required)" << endl; // TO DO: add depth estimation
                cout << " -k : k-mer length [17]" << endl;
                cout << " -a : use fixed alignment threshold [50]" << endl;
                cout << " -x : alignment x-drop factor [7]" << endl;
                cout << " -e : error rate [auto estimated from fastq]" << endl;
                cout << " -m : total RAM of the system in MB [auto estimated if possible or 8,000 if not]" << endl;
                cout << " -z : skip the pairwise alignment [false]" << endl;
                cout << " -w : relaxMargin parameter for alignment on edges [300]" << endl;
                cout << " -c : alignment score deviation from the mean [0.1]" << endl;
                cout << " -n : filter out alignment on edge [false]" << endl;
                cout << " -r : kmerRift: bases separating two k-mers used as seeds for a read [1,000]" << endl;
                cout << " -u : HOPC: use HOPC representation for kmers [false]" << endl; // TODO: pick a better letter because u doesn't make sense
                cout << " -p : output in PAF format [false]\n" << endl;

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
    int lower, upper; // reliable range lower and upper bound
    double ratioPhi;
    Kmer::set_k(kmer_len);
    size_t upperlimit = 10000000; // in bytes
    Kmers kmervect;
    vector<string> seqs;
    vector<string> quals;
    vector<string> nametags;
    readVector_ reads;
    Kmers kmersfromreads;
    vector<tuple<size_t,size_t,size_t>> occurrences;
    vector<tuple<size_t,size_t,size_t>> transtuples;

    // TODO: remove this test:
    cout << "\n\n\nTesting Kmer class" << endl;
    testKmer();
    cout << "\n\n\n" << endl;
    
    
    // 
    // File and setting used
    //
#ifdef PRINT
#ifdef JELLYFISH
    cout << "K-mer counting: Jellyfish" << endl;
    cout << "Input k-mer file: " << kmer_file << endl;
#endif
    cout << "K-mer counting: BELLA" << endl;
    cout << "Output filename: " << out_file << endl;
    cout << "K-mer length: " << kmer_len << endl;
    cout << "X-drop: " << xdrop << endl;
    cout << "Depth: " << depth << "X" << endl;
    if(b_parameters.skipAlignment)
        cout << "Compute alignment: false" << endl;
    else cout << "Compute alignment: true" << endl;
    if(!b_parameters.allKmer)
        cout << "Seeding: two-kmer" << endl;
    else cout << "Seeding: all-kmer" << endl;
#endif

    //
    // Kmer file parsing, error estimation, reliable bounds computation, and k-mer dictionary creation
    //

    dictionary_t countsreliable;
#ifdef JELLYFISH
    // Reliable bounds computation for Jellyfish using default error rate
    double all = omp_get_wtime();
    lower = computeLower(depth,erate,kmer_len);
    upper = computeUpper(depth,erate,kmer_len);
    cout << "Error rate is " << erate << endl;
    cout << "Reliable lower bound: " << lower << endl;
    cout << "Reliable upper bound: " << upper << endl;
    JellyFishCount(kmer_file, countsreliable, lower, upper);

#else
    // Error estimation and reliabe bounds computation within denovo counting
    cout << "\nRunning with up to " << MAXTHREADS << " threads" << endl;
    double all = omp_get_wtime();
    DeNovoCount(allfiles, countsreliable, lower, upper, kmer_len, depth, erate, upperlimit, b_parameters);

#ifdef PRINT
    cout << "Error rate estimate is " << erate << endl;
    cout << "Reliable lower bound: " << lower << endl;
    cout << "Reliable upper bound: " << upper << endl;
if(b_parameters.adapThr)
{
    ratioPhi = adaptiveSlope(erate);
    cout << "Deviation from expected alignment score: " << b_parameters.deltaChernoff << endl;
    cout << "Constant of adaptive threshold: " << ratioPhi*(1-b_parameters.deltaChernoff)<< endl;
}
else cout << "Default alignment score threshold: " << b_parameters.defaultThr << endl;
if(b_parameters.alignEnd)
{
    cout << "Constraint: alignment on edge with a margin of " << b_parameters.relaxMargin << " bps" << endl;
}
#endif // PRINT
#endif // DENOVO COUNTING

    //
    // Fastq(s) parsing
    //

    double parsefastq = omp_get_wtime();
    size_t read_id = 0; // read_id needs to be global (not just per file)

#ifdef PRINT
    cout << "\nRunning with up to " << MAXTHREADS << " threads" << endl;
#endif

    vector < vector<tuple<int,int,int>> > alloccurrences(MAXTHREADS);   
    vector < vector<tuple<int,int,int>> > alltranstuples(MAXTHREADS);   
    vector < readVector_ > allreads(MAXTHREADS);

    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++)
    {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);

        size_t fillstatus = 1;
        while(fillstatus)
        { 
            fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
            size_t nreads = seqs.size();

            #pragma omp parallel for
            for(int i=0; i<nreads; i++) 
            {
                // remember that the last valid position is length()-1
                int len = seqs[i].length();

                readType_ temp;
                nametags[i].erase(nametags[i].begin());     // removing "@"
                temp.nametag = nametags[i];
                temp.seq = seqs[i];     // save reads for seeded alignment
                temp.readid = read_id+i;

                allreads[MYTHREAD].push_back(temp);
                        
                for(int j=0; j<=len-kmer_len; j++)  
                {
                    std::string kmerstrfromfastq = seqs[i].substr(j, kmer_len);
                    Kmer mykmer(kmerstrfromfastq.c_str(), kmerstrfromfastq.length());
                    bool found;
                    int idx; // kmer_id
                    Kmer lexsmall;
                    // remember to use only ::rep() when building kmerdict as well
                    // TODO: figure out what kmerdict is, because it doesn't actually exist
                    if (b_parameters.useHOPC) { // TODO: see if this if should be on the outside of a loop for efficiency
                      lexsmall = mykmer.getHOPC(); // ::getHOPC() includes ::rep()
                    } else {
                      lexsmall = mykmer.rep();
                    }
                    found = countsreliable.find(lexsmall, idx);

                    if(found)
                    {
                        alloccurrences[MYTHREAD].emplace_back(std::make_tuple(read_id+i,idx,j)); // vector<tuple<read_id,kmer_id,kmerpos>>
                        alltranstuples[MYTHREAD].emplace_back(std::make_tuple(idx,read_id+i,j)); // transtuples.push_back(col_id,row_id,kmerpos)
                    }
                }
            } // for(int i=0; i<nreads; i++)
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

#ifdef PRINT
    cout << "Fastq(s) parsing fastq took: " << omp_get_wtime()-parsefastq << "s" << endl;
    cout << "Total number of reads: "<< read_id << "\n"<< endl;
#endif

    //
    // Sparse matrices construction
    //

    double matcreat = omp_get_wtime();

    size_t nkmer = countsreliable.size();
    CSC<size_t,size_t> spmat(occurrences, read_id, nkmer, 
                            [] (size_t & p1, size_t & p2) 
                            {  
                                return p1;
                            });
    std::vector<tuple<size_t,size_t,size_t>>().swap(occurrences); // remove memory of occurences

    CSC<size_t,size_t> transpmat(transtuples, nkmer, read_id, 
                            [] (size_t & p1, size_t & p2) 
                            {  return p1;
                            });
    std::vector<tuple<size_t,size_t,size_t>>().swap(transtuples); // remove memory of transtuples

#ifdef PRINT
    cout << "Sparse matrix construction took: " << omp_get_wtime()-matcreat << "s\n" << endl;
#endif

    //
    // Overlap detection (sparse matrix multiplication) and seed-and-extend alignment
    //

    spmatPtr_ getvaluetype(make_shared<spmatType_>());
    HashSpGEMM(spmat, transpmat, 
            [] (size_t & pi, size_t & pj)                     // n-th k-mer positions on read i and on read j
            {   spmatPtr_ value(make_shared<spmatType_>());
                value->count = 1;
                value->pos.push_back(make_pair(pi, pj));
                return value;
            },
    [&kmer_len,&b_parameters] (spmatPtr_ & m1, spmatPtr_ & m2)
            {
                for(int i = 0; i < m1->pos.size(); ++i)
                {
                    int left  = m2->pos[i].first - kmer_len - b_parameters.kmerRift;
                    int right = m2->pos[i].first + kmer_len + b_parameters.kmerRift;
                    int newseed  = m1->pos[i].first;

                    if(!isinrift(newseed, left, right))       // seeds separated by <kmerRift> bases
                    {
                        left  = m2->pos[i].second - kmer_len - b_parameters.kmerRift;
                        right = m2->pos[i].second + kmer_len + b_parameters.kmerRift;
                        newseed  = m1->pos[i].second;

                        if(!isinrift(newseed, left, right))   // seeds separated by <kmerRift> bases
                            if(!b_parameters.allKmer)         // save at most two kmers as seeds
                            {
                                m2->count = m2->count+m1->count;
                                m2->pos.clear();  // free from previous positions and save only two pos

                                m2->pos.push_back(make_pair(m2->pos[i].first, m2->pos[i].second));
                                m2->pos.push_back(make_pair(m1->pos[i].first, m1->pos[i].second));

                                break;
                            }
                            else   // save all possible kmers as seeds
                            {
                                m2->count = m2->count+m1->count;
                                m2->pos.push_back(m1->pos[i]);
                            }
                    }
                }
                return m2;
            }, reads, getvaluetype, kmer_len, xdrop, out_file, b_parameters, ratioPhi); 

    cout << "Total running time: " << omp_get_wtime()-all << "s\n" << endl;
    return 0;
} 
