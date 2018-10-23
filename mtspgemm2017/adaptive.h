#include "CSC.h"
#include "alignment.h"
#include "global.h"
#include "common.h"
#include "../kmercode/hash_funcs.h"
#include "../kmercode/Kmer.hpp"
#include "../kmercode/Buffer.h"
#include "../kmercode/common.h"
#include "../kmercode/fq_reader.h"
#include "../kmercode/ParallelFASTQ.h"
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/seeds.h>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_statistics.h>

using namespace seqan;

typedef Seed<Simple>  TSeed;
typedef SeedSet<TSeed> TSeedSet;

#define PERCORECACHE (1024 * 1024)
//#define TIMESTEP
//#define PRINT
//#define RAM
//#define OSX
//#define LINUX
//#define THREADLIMIT
//#define MAX_NUM_THREAD 1

#ifdef OSX
#include <mach/mach.h>
#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
#endif

#ifdef LINUX
#include "sys/types.h"
#include "sys/sysinfo.h"
struct sysinfo info;
#endif

/**
 * @brief writeToFile writes a CSV containing
 * CSC indices of the output sparse matrix
 * @param myBatch
 * @param filename
 */
void writeToFile(std::stringstream & myBatch, std::string filename)
{  
    std::string myString = myBatch.str();  
    std::ofstream myfile;  
    myfile.open (filename, ios_base::app);  
    myfile << myString;  
    myfile.close();  
}

bool onedge(int colStart, int colEnd, int colLen, int rowStart, int rowEnd, int rowLen)
{
    int minLeft = min(colStart, rowStart);
    int minRight = min(colLen-colEnd, rowLen-rowEnd);
    int epsilon = 300;

    if(minLeft-epsilon <= 0)
        minLeft = 0;

    if(minRight-epsilon <= 0)
        minRight = 0;

    if((minLeft == 0 && minRight == 0))
        return true;
    else
        return false;
}

template <typename IT, typename NT, typename MultiplyOperation, typename AddOperation, typename FT>
void LocalSpGEMM(IT & start, IT & end, IT & ncols, const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, vector<IT> * RowIdsofC, vector<FT> * ValuesofC)
{
    int v = 0;
    for(int i = start; i<end; ++i) // for bcols of B (one block)
    {
        IT hsize = B.colptr[i+1]-B.colptr[i];
        HeapEntry<IT,FT> *mergeheap = new HeapEntry<IT,FT>[hsize];
        IT maxnnzc = 0;           // max number of nonzeros in C(:,i)

        IT k = 0;   // make initial heap

        for(IT j=B.colptr[i]; j < B.colptr[i+1]; ++j) // for all the nonzeros of the ith column
        {
            IT inner = B.rowids[j]; // get the row id of B (or column id of A)
            IT npins = A.colptr[inner+1] - A.colptr[inner]; // get the number of nonzeros in A's corresponding column
        
            if(npins > 0)
            {
                mergeheap[k].loc = 1;
                mergeheap[k].runr = j; // the pointer to B.rowid's is the run-rank
                mergeheap[k].value = multop(A.values[A.colptr[inner]],B.values[j]);
                mergeheap[k++].key = A.rowids[A.colptr[inner]]; // A's first rowid is the first key
                maxnnzc += npins;
            }
        }

        hsize = k; // if any of A's "significant" columns is empty, k will be less than hsize
        make_heap(mergeheap, mergeheap + hsize);        
        // reserve changes the capacity of the vector, so that future push_back's won't cause reallocation
        // but it does not change the size, you can still use v.size() to get the number of valid elements
 
        RowIdsofC[v].reserve(maxnnzc); 
        ValuesofC[v].reserve(maxnnzc);
        
        while(hsize > 0)
        {
            pop_heap(mergeheap, mergeheap + hsize);         // result is stored in mergeheap[hsize-1]
            HeapEntry<IT,FT> hentry = mergeheap[hsize-1];
            
            // Use short circuiting
            if( (!RowIdsofC[v].empty()) && RowIdsofC[v].back() == hentry.key)
            {
                ValuesofC[v].back() = addop(hentry.value, ValuesofC[v].back());
            }
            else
            {
                ValuesofC[v].push_back(hentry.value);
                RowIdsofC[v].push_back(hentry.key);
            }
      
            IT inner = B.rowids[hentry.runr];
            
            // If still unused nonzeros exists in A(:,colind), insert the next nonzero to the heap
            if( (A.colptr[inner + 1] - A.colptr[inner]) > hentry.loc)
            {
                IT index = A.colptr[inner] + hentry.loc;
                mergeheap[hsize-1].loc   = hentry.loc +1;
                mergeheap[hsize-1].runr  = hentry.runr;
                mergeheap[hsize-1].value = multop(A.values[index], B.values[hentry.runr]);
                mergeheap[hsize-1].key   = A.rowids[index];
                
                push_heap(mergeheap, mergeheap + hsize);
            }
            else
            {
                --hsize;
            }     
        }
        v++;
        delete [] mergeheap;
    }
}

/**
  * Sparse multithreaded GEMM.
  * Probably slower than HeapSpGEMM_gmalloc but likely to use less memory
 **/
template <typename IT, typename NT, typename FT, typename MultiplyOperation, typename AddOperation>
void HeapSpGEMM(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, readVector_ & read, 
    FT & getvaluetype, int kmer_len, int algnmnt_drop, int algnmnt_thr, char* filename, bool skip_algnmnt_krnl)
{   
    size_t upperlimit = 10000000; // in bytes
#ifdef RAM // number of cols depends on available RAM
    cout << "Cols subdivision based on available RAM\n" << endl;

#ifdef OSX // OSX-based memory consumption implementation 
    vm_size_t page_size;
    mach_port_t mach_port;
    mach_msg_type_number_t count;
    vm_statistics64_data_t vm_stats;
    uint64_t free_memory, used_memory;

    mach_port = mach_host_self();
    count = sizeof(vm_stats) / sizeof(natural_t);

    if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
                KERN_SUCCESS == host_statistics64(mach_port, HOST_VM_INFO,
                                                (host_info64_t)&vm_stats, &count))
    {
        free_memory = (int64_t)vm_stats.free_count * (int64_t)page_size;
        used_memory = ((int64_t)vm_stats.active_count +
                      (int64_t)vm_stats.inactive_count +
                      (int64_t)vm_stats.wire_count) * (int64_t)page_size;
        /* It is good to note that just because Mac OS X may show very little actual free memory at times that it may 
        not be a good indication of how much is ready to be used on short notice. */
    } 
#endif
    
#ifdef LINUX // LINUX-based memory consumption implementation 
    if(sysinfo(&info) != 0)
    {
        return false;
    }   
    unsigned long free_memory = info.freeram * info.mem_unit;
    free_memory += info.freeswap * info.mem_unit;
    free_memory += info.bufferram * info.mem_unit;
#endif

    intmax_t flops = 0; // total flops (multiplication) needed to generate C
#pragma omp parallel 
    {
        intmax_t tflops=0; //thread private flops
#pragma omp for
        for(int i=0; i < B.cols; ++i)        // for all columns of B
        {
            IT locmax = 0;
            for(IT j=B.colptr[i]; j < B.colptr[i+1]; ++j)    // For all the nonzeros of the ith column
            {
        //cout << B.colptr[i] << " B.colptr[i]" << endl;
                IT inner = B.rowids[j];             // get the row id of B (or column id of A)
                IT npins = A.colptr[inner+1] - A.colptr[inner]; // get the number of nonzeros in A's corresponding column
        //cout << npins << " npins" << endl;
                locmax += npins;
            }

            tflops += (intmax_t)locmax; 
        }
#pragma omp critical
        {
            flops += tflops;
        }
    }
    cout << "Flops: " << flops << endl; 

    uint64_t d = B.nnz/B.cols; // about d nnz each col
    uint64_t rsv = B.nnz*d;    // worst case
    uint64_t required_memory = (rsv)*sizeof(FT);
    //IT * rowids;
    //FT * values;
    
    // here numThreads actually represents the number of blocks based on available RAM
    int numThreads = required_memory/free_memory; 
    int colsPerBlock = B.cols/numThreads;                 // define number of columns for each blocks

    // multi thread variable definition 
    int * colStart = new int[numThreads+1];              // need one block more for remaining cols
    int * colEnd = new int[numThreads+1];                // need one block more for remaining cols
    int * numCols = new int[numThreads+1];
#ifdef PRINT
    cout << "numBlocks " << numThreads+1 << endl;
    cout << "colsPerBlock " << colsPerBlock << "\n" << endl;
#endif
    colStart[0] = 0;
    colEnd[0] = 0;
    numCols[0] = 0;

    int colsTrace = 0;
    for(int i = 0; i < numThreads+1; ++i)
    {
        colStart[i] = colsTrace;
        if(colsTrace+colsPerBlock <= B.cols)
        {
            colEnd[i] = colStart[i]+colsPerBlock;
            numCols[i] = colsPerBlock;

        } else 
        {
            colEnd[i] = B.cols;
            numCols[i] = B.cols-colsTrace;
        }
        colsTrace = colsTrace+numCols[i];
    }
#else // number of cols depends on number of threads
    cout << "Cols subdivision based on number of threads\n" << endl;

    intmax_t flops = 0; // total flops (multiplication) needed to generate C
#pragma omp parallel 
    {
        intmax_t tflops=0; //thread private flops
#pragma omp for
        for(int i=0; i < B.cols; ++i)        // for all columns of B
        {
            IT locmax = 0;
            for(IT j=B.colptr[i]; j < B.colptr[i+1]; ++j)    // For all the nonzeros of the ith column
            {
        //cout << B.colptr[i] << " B.colptr[i]" << endl;
                IT inner = B.rowids[j];             // get the row id of B (or column id of A)
                IT npins = A.colptr[inner+1] - A.colptr[inner]; // get the number of nonzeros in A's corresponding column
        //cout << npins << " npins" << endl;
                locmax += npins;
            }

            tflops += (intmax_t)locmax; 
        }
#pragma omp critical
        {
            flops += tflops;
        }
    }
    cout << "Flops: " << flops << endl; 

#ifdef THREADLIMIT
    omp_set_dynamic(0);                      // Explicitly disable dynamic teams
    omp_set_num_threads(MAX_NUM_THREAD);     // Use MAX_NUM_THREAD threads for all consecutive parallel regions
#endif

    int numThreads;
#pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }

    //IT * rowids;
    //FT * values;

    int colsPerBlock = B.cols/numThreads;                // define number of columns for each blocks

    // multithread variable definition
    int * colStart = new int[numThreads+1];              // need one block more for remaining cols
    int * colEnd = new int[numThreads+1];                // need one block more for remaining cols
    int * numCols = new int[numThreads+1];
#ifdef PRINT
    cout << "numThreads: " << numThreads << endl;
    cout << "colsPerThread: " << colsPerBlock << "\n" << endl;
#endif
    colStart[0] = 0;
    colEnd[0] = 0;
    numCols[0] = 0;

    int colsTrace = 0;
    for(int i = 0; i < numThreads+1; ++i)
    {
        colStart[i] = colsTrace;
        if(colsTrace+colsPerBlock <= B.cols)
        {
            colEnd[i] = colStart[i]+colsPerBlock;
            numCols[i] = colsPerBlock;

        } else 
        {
            colEnd[i] = B.cols;
            numCols[i] = B.cols-colsTrace;
        }
        colsTrace = colsTrace+numCols[i];
    }
    #endif

    //IT * colptr = new IT[B.cols+1]; 
    //colptr[0] = 0;

    shared_ptr<readVector_> globalInstance = make_shared<readVector_>(read); 
    
    map<string,int> readset; // GGGG: reads for inner spgemm
    vector<pair<int,int>> idxset; // GGGG: indexes the pairs we care about in the spgemm (could be useless)
    stringstream myBatch; // each thread saves its results in its provate stringstream variable
    seqAnResult longestExtensionScore;
    double ovlalign = omp_get_wtime(); // get overlap and alignment time

#pragma omp parallel for private(myBatch, longestExtensionScore) shared(colStart,colEnd,numCols,globalInstance)
    for(int b = 0; b < numThreads+1; ++b) 
    {
#ifdef TIMESTEP
        double ovl = omp_get_wtime();
#endif
        //spmatPtr_ getvaluetype2(make_shared<spmatType_>());
        vector<IT> * RowIdsofC = new vector<IT>[numCols[b]];    // row ids for each column of C (bunch of cols)
        vector<FT> * ValuesofC = new vector<FT>[numCols[b]];    // values for each column of C (bunch of cols)

        IT * colptr = new IT[numCols[b]+1];
        colptr[0] = 0;

        LocalSpGEMM(colStart[b], colEnd[b], numCols[b], A, B, multop, addop, RowIdsofC, ValuesofC);

#ifdef TIMESTEP
#pragma omp critical
        {
            cout << "#" << omp_get_thread_num()+1 << ", ovelap time: " << omp_get_wtime()-ovl << "s" << endl;
        }
#endif
        int k=0;
        for(int i=0; i<numCols[b]; ++i) // for all edge lists (do NOT parallelize)
        {
            colptr[i+1] = colptr[i] + RowIdsofC[k].size();
            ++k;
        }

        IT * rowids = new IT[colptr[numCols[b]]];
        FT * values = new FT[colptr[numCols[b]]];

        k=0;
        for(int i=0; i<numCols[b]; ++i) // combine step
        {
            copy(RowIdsofC[k].begin(), RowIdsofC[k].end(), rowids + colptr[i]);
            copy(ValuesofC[k].begin(), ValuesofC[k].end(), values + colptr[i]);
            ++k;
        }

        intmax_t tot = 0;
        intmax_t tnum = 0;
        bool result;

        delete [] RowIdsofC;
        delete [] ValuesofC;

        // Local Alignment before write on stringstream 
        for(int i=0; i<numCols[b]; ++i) 
        {
            for(int j=colptr[i]; j<colptr[i+1]; ++j) 
            {
                if(values[j]->count == 1)
                {
                    longestExtensionScore = seqanAlOne(globalInstance->at(rowids[j]).seq, globalInstance->at(i+colStart[b]).seq, 
                        globalInstance->at(rowids[j]).seq.length(), values[j]->pos[0], values[j]->pos[1], algnmnt_drop, kmer_len);

                    int diffCol = endPositionV(longestExtensionScore.seed)-beginPositionV(longestExtensionScore.seed);
                    int diffRow = endPositionH(longestExtensionScore.seed)-beginPositionH(longestExtensionScore.seed);
                    int minLeft = min(beginPositionV(longestExtensionScore.seed), beginPositionH(longestExtensionScore.seed));
                    int minRight = min(globalInstance->at(i+colStart[b]).seq.length()-endPositionV(longestExtensionScore.seed), globalInstance->at(rowids[j]).seq.length()-endPositionH(longestExtensionScore.seed));

                    int ov = minLeft+minRight+(diffCol+diffRow)/2;
                    //double adaptive_thr = 0.4450*(double)ov;    // scoring matrix 1,-1,-1
                    //double adaptive_thr = 1.1675*(double)ov;    // scoring matrix 2,-1,-1
		//	double adaptive_thr = 0.3319*(double)ov+596.4349; //vv
//double adaptive_thr = 0.383265617669481*(double)ov+477.781646654574; //ab  
//double adaptive_thr = 1.10014413143354*(double)ov+714.940575445321; //ab                    
double adaptive_thr = 0.236875*(double)ov; 
//double adaptive_thr = 0.11*(double)ov + 165.38; // real ecoli dataset x = 3
                    //double adaptive_thr = 0.33*(double)ov + 286.86; // real ecoli dataset x = 7
                    //double adaptive_thr = 0.38*(double)ov + 215.66; // real ecoli dataset x = 9
                    //double adaptive_thr = 0.40*(double)ov + 195.17; // real ecoli dataset x = 11
                    //double adaptive_thr = 0.43*(double)ov + 138.37; // real ecoli dataset x = 13
                    //double adaptive_thr = 0.42*(double)ov + 44.57; // real ecoli dataset x = 15

                    //double adaptive_thr = 0.19*(double)ov + 451.41; // synthetic paeruginosa dataset x = 3
                    //double adaptive_thr = 0.12*(double)ov + 1047.16; // synthetic paeruginosa dataset x = 7
                    //double adaptive_thr = 0.11*(double)ov + 1059.30; // synthetic paeruginosa dataset x = 9
                    //double adaptive_thr = 0.18*(double)ov + 912.33; // synthetic paeruginosa dataset x = 11
                    //double adaptive_thr = 0.11*(double)ov + 1071.75; // synthetic paeruginosa dataset x = 13
                    //double adaptive_thr = 0.21*(double)ov + 828.65; // synthetic paeruginosa dataset x = 15

                    if((double)longestExtensionScore.score > adaptive_thr)
                    {
                    myBatch << globalInstance->at(i+colStart[b]).nametag << '\t' << globalInstance->at(rowids[j]).nametag << '\t' << values[j]->count << '\t' << longestExtensionScore.score << '\t' << longestExtensionScore.strand << '\t' << beginPositionV(longestExtensionScore.seed) << '\t' << 
                        endPositionV(longestExtensionScore.seed) << '\t' << globalInstance->at(i+colStart[b]).seq.length() << '\t' << beginPositionH(longestExtensionScore.seed) << '\t' << endPositionH(longestExtensionScore.seed) <<
                            '\t' << globalInstance->at(rowids[j]).seq.length() << endl;
                    }
                } // if(values[j]->count == 1)
                else
                {

                    longestExtensionScore = seqanAlGen(globalInstance->at(rowids[j]).seq, globalInstance->at(i+colStart[b]).seq, 
                        globalInstance->at(rowids[j]).seq.length(), values[j]->pos[0], values[j]->pos[1], values[j]->pos[2], values[j]->pos[3], algnmnt_drop, kmer_len);

                    int diffCol = endPositionV(longestExtensionScore.seed)-beginPositionV(longestExtensionScore.seed);
                    int diffRow = endPositionH(longestExtensionScore.seed)-beginPositionH(longestExtensionScore.seed);
                    int minLeft = min(beginPositionV(longestExtensionScore.seed), beginPositionH(longestExtensionScore.seed));
                    int minRight = min(globalInstance->at(i+colStart[b]).seq.length()-endPositionV(longestExtensionScore.seed), globalInstance->at(rowids[j]).seq.length()-endPositionH(longestExtensionScore.seed));

                    int ov = minLeft+minRight+(diffCol+diffRow)/2;
                    //double adaptive_thr = 0.4450*(double)ov;    // scoring matrix 1,-1,-1
                    //double adaptive_thr = 1.1675*(double)ov;    // scoring matrix 2,-1,-1
//double adaptive_thr = 0.383265617669481*(double)ov+477.781646654574; //ab 
//double adaptive_thr = 1.10014413143354*(double)ov+714.940575445321; //ab
//double adaptive_thr =0.8567*(double)ov+586.8543;  
double adaptive_thr = 0.236875*(double)ov;   
                //double adaptive_thr = 0.11*(double)ov + 165.38; // real ecoli dataset x = 3
                    //double adaptive_thr = 0.33*(double)ov + 286.86; // real ecoli dataset x = 7
                    //double adaptive_thr = 0.38*(double)ov + 215.66; // real ecoli dataset x = 9
                    //double adaptive_thr = 0.40*(double)ov + 195.17; // real ecoli dataset x = 11
                    //double adaptive_thr = 0.43*(double)ov + 138.37; // real ecoli dataset x = 13
                    //double adaptive_thr = 0.42*(double)ov + 44.57; // real ecoli dataset x = 15

                    //double adaptive_thr = 0.19*(double)ov + 451.41; // synthetic paeruginosa dataset x = 3
                    //double adaptive_thr = 0.12*(double)ov + 1047.16; // synthetic paeruginosa dataset x = 7
                    //double adaptive_thr = 0.11*(double)ov + 1059.30; // synthetic paeruginosa dataset x = 9
                    //double adaptive_thr = 0.18*(double)ov + 912.33; // synthetic paeruginosa dataset x = 11
                    //double adaptive_thr = 0.11*(double)ov + 1071.75; // synthetic paeruginosa dataset x = 13
                    //double adaptive_thr = 0.21*(double)ov + 828.65; // synthetic paeruginosa dataset x = 15

                    if((double)longestExtensionScore.score > adaptive_thr)
                    {

                    myBatch << globalInstance->at(i+colStart[b]).nametag << '\t' << globalInstance->at(rowids[j]).nametag << '\t' << values[j]->count << '\t' << longestExtensionScore.score << '\t' << longestExtensionScore.strand << '\t' << beginPositionV(longestExtensionScore.seed) << '\t' << 
                        endPositionV(longestExtensionScore.seed) << '\t' << globalInstance->at(i+colStart[b]).seq.length() << '\t' << beginPositionH(longestExtensionScore.seed) << '\t' << endPositionH(longestExtensionScore.seed) <<
                            '\t' << globalInstance->at(rowids[j]).seq.length() << endl;
                    }

                }
            }// for(int j=colptr[i]; j<colptr[i+1]; ++j)
        } // for(int i=0; i<numCols[b]; ++i)

#ifdef TIMESTEP
#pragma omp critical
        {
            cout << "#" << omp_get_thread_num()+1 << ", alignment time: " << omp_get_wtime()-align << "s" << endl;
        }
#endif

        delete [] colptr;
        delete [] rowids;
        delete [] values;

#pragma omp critical
        {
            writeToFile(myBatch, filename);
            myBatch.str(std::string());
        }
    }//for(int b = 0; b < numThreads+1; ++b)

    delete [] colStart;
    delete [] colEnd;
    delete [] numCols; 
}

