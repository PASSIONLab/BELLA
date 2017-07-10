#include "CSC.h"
#include "global.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <mach/mach.h>
#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>

/* OSX-BASED BLOCK MULTIPLICATION */

#define PERCORECACHE (1024 * 1024)
#define KMER_LENGTH 17

template <typename IT, typename NT, typename FT, typename MultiplyOperation, typename AddOperation>
void LocalSpGEMM(IT & start, IT & end, IT & ncols, const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, vector<IT> * RowIdsofC, vector<FT> * ValuesofC)
{
    int i, v;
    //#pragma omp parallel for
    for(i = start, v=0; i<end, v<ncols; ++i, v++) // for bcols of B (one block)
    {
        IT hsize = B.colptr[i+1]-B.colptr[i];
        HeapEntry<IT,FT> *mergeheap = new HeapEntry<IT,FT>[hsize];
        IT maxnnzc = 0;      // max number of nonzeros in C(:,i)

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
        delete [] mergeheap;
    }

}

/**
  * Sparse multithreaded GEMM.
  * Probably slower than HeapSpGEMM_gmalloc but likely to use less memory
 **/
//template <typename IT, typename NT, typename FT, typename MultiplyOperation, typename AddOperation, typename AlphaOperation>
template <typename IT, typename NT, typename FT, typename MultiplyOperation, typename AddOperation>
void HeapSpGEMM(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, CSC<IT,FT> & C)
{
    IT start, end, ncols;
    
    /* OSX-based memory consumption implementation */

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
        
        cout << "Free memory: " << free_memory/PERCORECACHE << "\nUsed memory: " << used_memory/PERCORECACHE << endl;
        /* It is good to note that just because Mac OS X may show very little actual free memory at times that it may 
        not be a good indication of how much is ready to be used on short notice. */

    } 

    uint64_t d = B.nnz/B.cols; // about d nnz each col
    uint64_t rsv = B.nnz*d;    // worst case
    uint64_t requiredmemory = (rsv)*sizeof(size_t);

    cout << "Required memory: " << requiredmemory/PERCORECACHE << endl;
    
    if(requiredmemory > free_memory)
    {
        cout << "*** BLOCK MULTIPLICATION ***" << endl; 
        IT blocks = requiredmemory/free_memory;
        IT bcols = B.cols/blocks;
        
        cout << "blocks = " << blocks << endl;
        cout << "bcols = " << bcols << endl;    // define number of columns for each blocks

        C.rows = A.rows;
        C.cols = B.cols;
        C.colptr = new IT[C.cols+1];
        C.colptr[0] = 0;

        IT pcols = 0;

        for(IT b = 0; b < blocks+1; b++) 
        { 
            start = pcols;
            if(pcols+bcols <= B.cols)
            {
                end = pcols+bcols;
                ncols = bcols;
            } else 
            {
                end = B.cols;
                ncols = B.cols-pcols;
            }
            
            vector<IT> * RowIdsofC = new vector<IT>[ncols];      // row ids for each column of C (bunch of cols)
            vector<FT> * ValuesofC = new vector<FT>[ncols];      // values for each column of C (bunch of cols)
        
            LocalSpGEMM(start, end, ncols, A, B, multop, addop, RowIdsofC, ValuesofC);
            
            int i, j;
            for(i=start, j=0; i<end, j<ncols; ++i, j++) // for all edge lists (do NOT parallelize)
            {
                C.colptr[i+1] = C.colptr[i] + RowIdsofC[j].size();
            }

            if(b == 0) // first iteration on blocks
            {
                C.nnz = C.colptr[end];
                C.rowids = new IT[C.nnz];
                C.values = new FT[C.nnz];

            } else // resize
            {
                IT nnznew = C.colptr[end];

                IT *temprow = new IT[nnznew];
                copy(C.rowids, C.rowids+C.nnz, temprow);
                delete C.rowids;

                C.rowids = new IT[nnznew];
                copy(temprow, temprow+nnznew, C.rowids);
                delete temprow;

                FT *tempval = new IT[nnznew];
                copy(C.values, C.values+C.nnz, tempval);
                delete C.values;

                C.values = new FT[nnznew];
                copy(tempval, tempval+nnznew, C.values);
                delete tempval;

                C.nnz = nnznew;
            }
         
            int m, l;
            //#pragma omp parallel for
            for(m=start, l=0; m<end, l<ncols; ++m, l++) // combine step
            {
                // to do: local alignment as alpha operation
                // transform(&ValuesofC[colStart[i]], &ValuesofC[colEnd[i]], &ValuesofC[colStart[i]], alphaop);
                copy(RowIdsofC[l].begin(), RowIdsofC[l].end(), C.rowids + C.colptr[m]);
                copy(ValuesofC[l].begin(), ValuesofC[l].end(), C.values + C.colptr[m]);
            }

            pcols = pcols+ncols; // update counter
            delete [] RowIdsofC;
            delete [] ValuesofC;
        }
    } else
    {
        cout << "*** STANDARD MULTIPLICATION ***" << endl;
        start = 0;
        end = B.cols;
        ncols = B.cols;
        vector<IT> * RowIdsofC = new vector<IT>[B.cols];      // row ids for each column of C
        vector<FT> * ValuesofC = new vector<FT>[B.cols];      // values for each column of C

        LocalSpGEMM(start, end, ncols, A, B, multop, addop, RowIdsofC, ValuesofC);

        if(C.isEmpty())
        {
            C.rows = A.rows;
            C.cols = B.cols;
            C.colptr = new IT[C.cols+1];
            C.colptr[0] = 0;
        
            for(int i=0; i < C.cols; ++i)        // for all edge lists (do NOT parallelize)
            {
                C.colptr[i+1] = C.colptr[i] + RowIdsofC[i].size();
            }

            C.nnz = C.colptr[C.cols];
            C.rowids = new IT[C.nnz];
            C.values = new FT[C.nnz];
        
            #pragma omp parallel for
            for(int i=0; i<C.cols; ++i)         // combine step
            {
                // to do: local alignment as alpha operation
                // transform(&ValuesofC[colStart[i]], &ValuesofC[colEnd[i]], &ValuesofC[colStart[i]], alphaop);
                copy(RowIdsofC[i].begin(), RowIdsofC[i].end(), C.rowids + C.colptr[i]);
                copy(ValuesofC[i].begin(), ValuesofC[i].end(), C.values + C.colptr[i]);
            }

            delete [] RowIdsofC;
            delete [] ValuesofC;
        }
    } 
}

/**
 ** with global memory allocation
 **/
template <typename IT, typename NT, typename FT, typename MultiplyOperation, typename AddOperation>
void HeapSpGEMM_gmalloc(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, CSC<IT,FT> & C )
{
    
    int numThreads;
#pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }
    
    // *************** Creating global space to store result, used by all threads *********************
    IT* maxnnzc = new IT[B.cols]; // maximum number of nnz in each column of C
    IT flops = 0; // total flops (multiplication) needed to generate C
#pragma omp parallel
    {
        IT tflops=0; //thread private flops
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

            maxnnzc[i] = locmax;
            tflops += locmax; 
        }
#pragma omp critical
        {
            flops += tflops;
        }
    } 
    IT flopsPerThread = flops/numThreads; // amount of work that will be assigned to each thread
    IT colPerThread [numThreads + 1]; // thread i will process columns from colPerThread[i] to colPerThread[i+1]-1
    
    IT* colStart = new IT[B.cols]; //start index in the global array for storing ith column of C
    IT* colEnd = new IT[B.cols]; //end index in the global array for storing ith column of C
    colStart[0] = 0;
    colEnd[0] = 0;
    
    int curThread = 0;
    colPerThread[curThread++] = 0;
    IT nextflops = flopsPerThread;
    
    // TODO: the following prefix sum can be parallelized, e.g., see
    // http://stackoverflow.com/questions/21352816/prefix-sums-taking-too-long-openmp
    // not a dominating term at this moment
    for(int i=0; i < (B.cols -1); ++i)
    {
        colStart[i+1] = colStart[i] + maxnnzc[i];
        colEnd[i+1] = colStart[i+1];
        if(nextflops < colStart[i+1])
        {
            colPerThread[curThread++] = i+1;
            nextflops += flopsPerThread;
        }
    }
    while(curThread < numThreads)
        colPerThread[curThread++] = B.cols;
    colPerThread[numThreads] = B.cols;

    
    
    IT size = colEnd[B.cols-1] + maxnnzc[B.cols-1];
    IT * RowIdsofC = new IT[size];
    FT * ValuesofC = new FT[size];
    delete [] maxnnzc;
    // ************************ End Creating global space *************************************
    
    
    // *************** Creating global heap space to be used by all threads *********************

    IT threadHeapSize[numThreads];
#pragma omp parallel
    {
        int thisThread = omp_get_thread_num();
        IT localmax = 0; // change -1 to 0 as our IT is size_t
        for(int i=colPerThread[thisThread]; i < colPerThread[thisThread+1]; ++i)
        {
            IT colnnz = B.colptr[i+1]-B.colptr[i];
            if(colnnz > localmax) 
        localmax = colnnz;
        }
        threadHeapSize[thisThread] = localmax;
    }
    
    IT threadHeapStart[numThreads+1];
    threadHeapStart[0] = 0;
    for(int i=0; i<numThreads; i++) {
        threadHeapStart[i+1] = threadHeapStart[i] + threadHeapSize[i];
    }
    
    HeapEntry<IT,FT> * globalheap = new HeapEntry<IT,FT>[threadHeapStart[numThreads]];
    
    // ************************ End Creating global heap space *************************************
    
#pragma omp parallel
    {
        int thisThread = omp_get_thread_num();
        HeapEntry<IT,FT> * mergeheap = globalheap + threadHeapStart[thisThread]; // thread private heap space
        
        
        for(int i=colPerThread[thisThread]; i < colPerThread[thisThread+1]; ++i)
        {
            IT k = 0;   // Make initial heap
            for(IT j=B.colptr[i]; j < B.colptr[i+1]; ++j)    // For all the nonzeros of the ith column
            {
                IT inner = B.rowids[j];             // get the row id of B (or column id of A)
                IT npins = A.colptr[inner+1] - A.colptr[inner]; // get the number of nonzeros in A's corresponding column
                
                if(npins > 0)
                {
                    mergeheap[k].loc = 1;
                    mergeheap[k].runr = j;              // the pointer to B.rowid's is the run-rank
                    mergeheap[k].value = multop(A.values[A.colptr[inner]],B.values[j]);
                    mergeheap[k++].key = A.rowids[A.colptr[inner]]; // A's first rowid is the first key
                }
            }
            IT hsize = k;      // if any of A's "significant" columns is empty, k will be less than hsize
            make_heap(mergeheap, mergeheap + hsize);
            
            // reserve changes the capacity of the vector, so that future push_back's won't cause reallocation
            // but it does not change the size, you can still use v.size() to get the number of valid elements
            // RowIdsofC[i].reserve(maxnnzc);
            // ValuesofC[i].reserve(maxnnzc);
            
            while(hsize > 0)
            {
                pop_heap(mergeheap, mergeheap + hsize);         // result is stored in mergeheap[hsize-1]
                HeapEntry<IT,FT> hentry = mergeheap[hsize-1];
                
                // Use short circuiting
                if( (colEnd[i] > colStart[i]) && RowIdsofC[colEnd[i]-1] == hentry.key)
                {
                    ValuesofC[colEnd[i]-1] = addop(hentry.value, ValuesofC[colEnd[i]-1]);
                }
                else
                {
                    
                    ValuesofC[colEnd[i]]= hentry.value ;
                    RowIdsofC[colEnd[i]]= hentry.key ;
                    colEnd[i] ++;
                }
                
                IT inner = B.rowids[hentry.runr];
                
                // If still unused nonzeros exists in A(:,colind), insert the next nonzero to the heap
                if( (A.colptr[inner + 1] - A.colptr[inner]) > hentry.loc)
                {
                    IT index = A.colptr[inner] + hentry.loc;
                    mergeheap[hsize-1].loc = hentry.loc +1;
                    mergeheap[hsize-1].runr = hentry.runr;
                    mergeheap[hsize-1].value = multop(A.values[index], B.values[hentry.runr]);
                    mergeheap[hsize-1].key = A.rowids[index];
                    
                    push_heap(mergeheap, mergeheap + hsize);
                }
                else
                {
                    --hsize;
                }
            }
        }
    }
    
    delete [] globalheap;
    
    if(C.isEmpty())
    {
        C.rows = A.rows;
        C.cols = B.cols;
        C.colptr = new IT[C.cols+1];
        C.colptr[0] = 0;
        
        for(int i=0; i < C.cols; ++i)        // for all edge lists (do NOT parallelize)
        {
            C.colptr[i+1] = C.colptr[i] + colEnd[i]-colStart[i];
        }

        C.nnz = C.colptr[C.cols];
        C.rowids = new IT[C.nnz];
        C.values = new FT[C.nnz];
        
#pragma omp parallel for
        for(int i=0; i< C.cols; ++i)         // combine step
        {
            copy(&RowIdsofC[colStart[i]], &RowIdsofC[colEnd[i]], C.rowids + C.colptr[i]);
            copy(&ValuesofC[colStart[i]], &ValuesofC[colEnd[i]], C.values + C.colptr[i]);
        }
  }
    delete [] RowIdsofC;
    delete [] ValuesofC;
    delete [] colEnd;
    delete [] colStart;
}