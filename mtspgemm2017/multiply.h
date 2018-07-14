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

using namespace seqan;

typedef Seed<Simple>  TSeed;
typedef SeedSet<TSeed> TSeedSet;

#define PERCORECACHE (1024 * 1024)
#define SEQAN
//#define TIMESTEP
//#define PRINT
//#define RAM
//#define OSX
//#define LINUX
//#define THREADLIMIT
//#define MAX_NUM_THREAD 4
//#define ALLKMER

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
    IT * rowids;
    FT * values;
    
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

    IT * rowids;
    FT * values;
    
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

    shared_ptr<readVector_> globalInstance = make_shared<readVector_>(read); // shared pointer to access readVector_ in omp safe way
    stringstream myBatch;                                                    // each thread saves its results in its provate stringstream variable
    seqAnResult longestExtensionScore;
    double ovlalign = omp_get_wtime(); // get overlap and alignment time
    intmax_t novl = 0; // debug counting overlaps

#pragma omp parallel for private(myBatch, longestExtensionScore) shared(colStart,colEnd,numCols,globalInstance)
    for(int b = 0; b < numThreads+1; ++b) 
    { 
#ifdef TIMESTEP
        double ovl = omp_get_wtime();
#endif
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

        intmax_t tovl = 0; // debug counting overlaps per thread

        delete [] RowIdsofC;
        delete [] ValuesofC;

        // Local Alignment before write on stringstream 
        for(int i=0; i<numCols[b]; ++i) 
        {
            for(int j=colptr[i]; j<colptr[i+1]; ++j) 
            {
                ++tovl; // DEBUG
                int reduced_kmer_len = 7;

                //
                // De novo kmer counting
                //

                dictionary_t countskmerpair;
                DeNovoCountPair(globalInstance->at(rowids[j]).seq, globalInstance->at(i+colStart[b]).seq, countskmerpair, reduced_kmer_len, upperlimit);
    
                vector<tuple<int,int,int>> alloccurrences;
                vector<tuple<int,int,int>> alltranstuples;

                int len_row = globalInstance->at(rowids[j]).seq.length();
                int len_col = globalInstance->at(i+colStart[b]).seq.length();

                for(int k=0; k<=len_row-reduced_kmer_len; k++)  
                {
                    std::string kmerstrfromrow = globalInstance->at(rowids[j]).seq.substr(k, reduced_kmer_len);
                    Kmer mykmer(kmerstrfromrow.c_str());
                    // remember to use only ::rep() when building kmerdict as well
                    Kmer lexsmall = mykmer.rep();

                    int idx; // kmer_id
                    auto found = countskmerpair.find(lexsmall,idx);
                    if(found)
                    {
                        alloccurrences.emplace_back(std::make_tuple(j,idx,k)); // vector<tuple<read_id,kmer_id,kmerpos>>
                        alltranstuples.emplace_back(std::make_tuple(idx,j,k)); // transtuples.push_back(col_id,row_id,kmerpos)
                    }
                }

                for(int k=0; k<=len_col-reduced_kmer_len; k++)  
                {
                    std::string kmerstrfromcol = globalInstance->at(i+colStart[b]).seq.substr(k, reduced_kmer_len);
                    Kmer mykmer(kmerstrfromcol.c_str());
                    // remember to use only ::rep() when building kmerdict as well
                    Kmer lexsmall = mykmer.rep();
    
                    int idx; // kmer_id
                    auto found = countskmerpair.find(lexsmall,idx);
                    if(found)
                    {
                        alloccurrences.emplace_back(std::make_tuple(i,idx,k)); // vector<tuple<read_id,kmer_id,kmerpos>>
                        alltranstuples.emplace_back(std::make_tuple(idx,i,k)); // transtuples.push_back(col_id,row_id,kmerpos)
                    }
                }

                vector<tuple<int,int,int>> occurrences;
                vector<tuple<int,int,int>> transtuples;
    
                size_t tuplecount = alloccurrences.size();
    
                occurrences.resize(tuplecount);
                transtuples.resize(tuplecount);

                //
                // Sparse matrices construction
                //

                CSC<int, int> spmat(occurrences, 2, countskmerpair.size(), 
                    [] (int & p1, int & p2) 
                    { return p1; });
                std::cout << "spmat created with " << spmat.nnz << " nonzeros" << endl;
                std::vector<tuple<int,int,int>>().swap(occurrences);    // remove memory of occurences

                CSC<int, int> transpmat(transtuples, countskmerpair.size(), 2, 
                    [] (int & p1, int & p2) 
                    { return p1; });

                std::vector<tuple<int,int,int>>().swap(transtuples); // remove memory of transtuples

                spmat.Sorted();
                transpmat.Sorted();

                //
                // Sparse matrix multiplication
                //

                spmatPtr_ getvaluetype(make_shared<spmatType_>());
                PairHeapSpGEMM(spmat, transpmat, 
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
                }, read, getvaluetype, reduced_kmer_len);
            }
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
            novl += tovl;
        }

#pragma omp critical
        {
            writeToFile(myBatch, filename);
            myBatch.str(std::string());
        }
    }
#ifdef PRINT
    cout << "nOverlap: " << novl << endl;
    cout << "nAlignment: " << naln << endl;
    cout << "Ovelap detection and Alignment time: " << omp_get_wtime()-ovlalign << "s" << endl;
#endif

    delete [] colStart;
    delete [] colEnd;
    delete [] numCols; 
}

template <typename IT, typename NT, typename MultiplyOperation, typename AddOperation, typename FT>
void PairLocalSpGEMM(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, vector<IT> * RowIdsofC, vector<FT> * ValuesofC)
{
    for(int i = 0; i<B.cols; ++i) // for all columns of B
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
 
        RowIdsofC[i].reserve(maxnnzc); 
        ValuesofC[i].reserve(maxnnzc);
        
        while(hsize > 0)
        {
            pop_heap(mergeheap, mergeheap + hsize);         // result is stored in mergeheap[hsize-1]
            HeapEntry<IT,FT> hentry = mergeheap[hsize-1];
            
            // Use short circuiting
            if( (!RowIdsofC[i].empty()) && RowIdsofC[i].back() == hentry.key)
            {
                ValuesofC[i].back() = addop(hentry.value, ValuesofC[i].back());
            }
            else
            {
                ValuesofC[i].push_back(hentry.value);
                RowIdsofC[i].push_back(hentry.key);
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
        delete [] mergeheap;
    }
} 

/**
  * Sparse multithreaded GEMM.
  * Probably slower than HeapSpGEMM_gmalloc but likely to use less memory
 **/
template <typename IT, typename NT, typename FT, typename MultiplyOperation, typename AddOperation>
void PairHeapSpGEMM(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, readVector_ & read, FT & getvaluetype, int kmer_len)
{
    IT * colptr = new IT[B.cols+1];
    colptr[0] = 0;

    vector<IT> * RowIdsofC = new vector<IT>[B.cols];      // row ids for each column of C
    vector<FT> * ValuesofC = new vector<FT>[B.cols];      // values for each column of C

    PairLocalSpGEMM(A, B, multop, addop, RowIdsofC, ValuesofC);

    for(int i=0; i < B.cols; ++i)        // for all edge lists (do NOT parallelize)
    {
        colptr[i+1] = colptr[i] + RowIdsofC[i].size();
    }

    IT * rowids = new IT[B.cols+1];
    FT * values = new FT[B.cols+1];

    for(int i=0; i < B.cols; ++i)  // combine step
    {
        copy(RowIdsofC[i].begin(), RowIdsofC[i].end(), rowids + colptr[i]);
        copy(ValuesofC[i].begin(), ValuesofC[i].end(), values + colptr[i]);
    }

    delete [] RowIdsofC;
    delete [] ValuesofC;

    delete rowids;
    delete values;
}

