#include "CSC.h"
#include <omp.h>
//#include <tbb/scalable_allocator.h>

#define PERCORECACHE (1024 * 1024)




template <typename IT, typename NT, typename MultiplyOperation, typename AddOperation>
void HeapSpGEMM(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, vector<IT> * RowIdsofC, vector<NT> * ValuesofC)
{
    #pragma omp parallel for
    for(int i=0; i < B.cols; ++i)        // for all columns of B
    {
        IT hsize = B.colptr[i+1]-B.colptr[i];
        HeapEntry<IT,NT> * mergeheap = new HeapEntry<IT,NT>[hsize];
        IT maxnnzc = 0;      // max number of nonzeros in C(:,i)
        
        IT k = 0;   // Make initial heap
        for(IT j=B.colptr[i]; j < B.colptr[i+1]; ++j)    // For all the nonzeros of the ith column
        {
            IT inner = B.rowids[j];				// get the row id of B (or column id of A)
            IT npins = A.colptr[inner+1] - A.colptr[inner];	// get the number of nonzeros in A's corresponding column
            
            if(npins > 0)
            {
                mergeheap[k].loc = 1;
                mergeheap[k].runr = j;    			// the pointer to B.rowid's is the run-rank
                mergeheap[k].value = multop(A.values[A.colptr[inner]],B.values[j]);
                mergeheap[k++].key = A.rowids[A.colptr[inner]];	// A's first rowid is the first key
                maxnnzc += npins;
            }
        }
        hsize = k;      // if any of A's "significant" columns is empty, k will be less than hsize
        make_heap(mergeheap, mergeheap + hsize);
        
        // reserve changes the capacity of the vector, so that future push_back's won't cause reallocation
        // but it does not change the size, you can still use v.size() to get the number of valid elements
        RowIdsofC[i].reserve(maxnnzc);
        ValuesofC[i].reserve(maxnnzc);
        
        while(hsize > 0)
        {
            pop_heap(mergeheap, mergeheap + hsize);         // result is stored in mergeheap[hsize-1]
            HeapEntry<IT,NT> hentry = mergeheap[hsize-1];
            
            // Use short circuiting
            if( (!RowIdsofC[i].empty()) && RowIdsofC[i].back() == hentry.key)
            {
                ValuesofC[i].back() = addop(hentry.value, ValuesofC[i].back());
            }
            else
            {
                ValuesofC[i].push_back( hentry.value );
                RowIdsofC[i].push_back( hentry.key );
            }
            
            IT inner = B.rowids[hentry.runr];
            
            // If still unused nonzeros exists in A(:,colind), insert the next nonzero to the heap
            if( (A.colptr[inner + 1] - A.colptr[inner]) > hentry.loc)
            {
                IT index = A.colptr[inner] + hentry.loc;
                mergeheap[hsize-1].loc	 = hentry.loc +1;
                mergeheap[hsize-1].runr  = hentry.runr;
                mergeheap[hsize-1].value = multop(A.values[index], B.values[hentry.runr]);
                mergeheap[hsize-1].key	 = A.rowids[index];
                
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
  * Functionally more general than C = \alpha A*B + \beta C
  * AlphaOperation generalizes scaling of A*B
  *     alphaop = bind2nd(std::plus<NT>(), alpha) for constant alpha, then we have the BLAS usage
  *     alphaop = bind2nd(std::less<NT>(), threshold)?_1:0 for constant threshold, then it works like a drop threshold
 **/
template <typename IT, typename NT, typename MultiplyOperation, typename AddOperation, typename AlphaOperation>
void LocalSpGEMM(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, AlphaOperation alphaop, CSC<IT,NT> & C )
{
    vector<IT> * RowIdsofC = new vector<IT>[B.cols];      // row ids for each column of C
    vector<NT> * ValuesofC = new vector<NT>[B.cols];      // values for each column of C

    
    HeapSpGEMM(A, B, multop, addop, RowIdsofC, ValuesofC);
    
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
        C.values = new NT[C.nnz];
	
        #pragma omp parallel for
        for(int i=0; i< C.cols; ++i)         // combine step
        {
            transform(ValuesofC[i].begin(), ValuesofC[i].end(), ValuesofC[i].begin(), alphaop);
            copy(RowIdsofC[i].begin(), RowIdsofC[i].end(), C.rowids + C.colptr[i]);
            copy(ValuesofC[i].begin(), ValuesofC[i].end(), C.values + C.colptr[i]);
        }
        delete [] RowIdsofC;
        delete [] ValuesofC;
    }
}





/**
 ** with global memory allocation
 **/
template <typename IT, typename NT, typename MultiplyOperation, typename AddOperation, typename AlphaOperation>
void HeapSpGEMM_gmalloc(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, AlphaOperation alphaop, CSC<IT,NT> & C )
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
                IT inner = B.rowids[j];				// get the row id of B (or column id of A)
                IT npins = A.colptr[inner+1] - A.colptr[inner];	// get the number of nonzeros in A's corresponding column
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
    NT * ValuesofC = new NT[size];
    delete [] maxnnzc;
    // ************************ End Creating global space *************************************
    
    
    // *************** Creating global heap space to be used by all threads *********************

    IT threadHeapSize[numThreads];
#pragma omp parallel
    {
        int thisThread = omp_get_thread_num();
        IT localmax = -1;
        for(int i=colPerThread[thisThread]; i < colPerThread[thisThread+1]; ++i)
        {
            IT colnnz = B.colptr[i+1]-B.colptr[i];
            if(colnnz > localmax) localmax = colnnz;
        }
        threadHeapSize[thisThread] = localmax;
    }
    
    IT threadHeapStart[numThreads+1];
    threadHeapStart[0] = 0;
    for(int i=0; i<numThreads; i++)
        threadHeapStart[i+1] = threadHeapStart[i] + threadHeapSize[i];
    
    HeapEntry<IT,NT> * globalheap = new HeapEntry<IT,NT>[threadHeapStart[numThreads]];
    
    // ************************ End Creating global heap space *************************************
    

    
    
#pragma omp parallel
    {
        int thisThread = omp_get_thread_num();
        HeapEntry<IT,NT> * mergeheap = globalheap + threadHeapStart[thisThread]; // thread private heap space
        
        
        for(int i=colPerThread[thisThread]; i < colPerThread[thisThread+1]; ++i)
        {
            IT k = 0;   // Make initial heap
            for(IT j=B.colptr[i]; j < B.colptr[i+1]; ++j)    // For all the nonzeros of the ith column
            {
                IT inner = B.rowids[j];				// get the row id of B (or column id of A)
                IT npins = A.colptr[inner+1] - A.colptr[inner];	// get the number of nonzeros in A's corresponding column
                
                if(npins > 0)
                {
                    mergeheap[k].loc = 1;
                    mergeheap[k].runr = j;    			// the pointer to B.rowid's is the run-rank
                    mergeheap[k].value = multop(A.values[A.colptr[inner]],B.values[j]);
                    mergeheap[k++].key = A.rowids[A.colptr[inner]];	// A's first rowid is the first key
                }
            }
            IT hsize = k;      // if any of A's "significant" columns is empty, k will be less than hsize
            make_heap(mergeheap, mergeheap + hsize);
            
            
            // reserve changes the capacity of the vector, so that future push_back's won't cause reallocation
            // but it does not change the size, you can still use v.size() to get the number of valid elements
            //RowIdsofC[i].reserve(maxnnzc);
            //ValuesofC[i].reserve(maxnnzc);
            
            while(hsize > 0)
            {
                pop_heap(mergeheap, mergeheap + hsize);         // result is stored in mergeheap[hsize-1]
                HeapEntry<IT,NT> hentry = mergeheap[hsize-1];
                
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
                    mergeheap[hsize-1].loc	 = hentry.loc +1;
                    mergeheap[hsize-1].runr  = hentry.runr;
                    mergeheap[hsize-1].value = multop(A.values[index], B.values[hentry.runr]);
                    mergeheap[hsize-1].key	 = A.rowids[index];
                    
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
        C.values = new NT[C.nnz];
        
#pragma omp parallel for
        for(int i=0; i< C.cols; ++i)         // combine step
        {
            transform(&ValuesofC[colStart[i]], &ValuesofC[colEnd[i]], &ValuesofC[colStart[i]], alphaop);
            copy(&RowIdsofC[colStart[i]], &RowIdsofC[colEnd[i]], C.rowids + C.colptr[i]);
            copy(&ValuesofC[colStart[i]], &ValuesofC[colEnd[i]], C.values + C.colptr[i]);
        }
        
    }
    
    delete [] RowIdsofC;
    delete [] ValuesofC;
    delete [] colEnd;
    delete [] colStart;
}

