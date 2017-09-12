#include "CSC.h"
#include "global.h"
#include "metric.h"
//#include "../edlib/edlib/include/edlib.h"
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
#include <seqan/seeds.h>

#define PERCORECACHE (1024 * 1024)
#define KMER_LENGTH 17
#define _OSX
#define TETA 3
#define OMEGA 2*0.15*300
//#define _EDLIB
//#define _MULTPTR

#ifdef _OSX
#include <mach/mach.h>
#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
#endif

#ifdef _LINUX
#include "sys/types.h"
#include "sys/sysinfo.h"
struct sysinfo info;
#endif

#ifdef _EDLIB
#define pCORR .7225
#define pWRON .2775
#define minLEN 17
/* Compute potential overlap */
void getOverlaplen(std::pair<int, int> & values, std::string & row, std::string col, int & len, int & left) 
{
    /* Left margin passed as argument */
    int right;           /* Define the right margin */
    int overlaplength;   /* Expected overlap region length */
    int temp1, temp2;
    
    if(values.first <= values.second)
        left = values.first;
    else left = values.second;

    temp1 = row.length() - values.first;   /* 1st value referred to A's row */
    temp2 = col.length() - values.second;  /* 2nd value referred to B's col */

    if(temp1 <= temp2)
        right = temp1;
    else right = temp2;

    len = left+right; /* Estimated overlap */
}

/* EDLIB Local Alignment */
bool edlibOp(std::pair<int, int> & values, std::string & row, std::string & col, int & len, int & ed) 
{
    bool align;
    int left = 0;
    /* Compute the overlap potential length and obtain substring of 
    overlaplength from the considered pair and pass them to alignop */
    getOverlaplen(values, row, col, len, left);
    /* Obtain row and col substr of overlaplength */
    std::string substrow = row.substr(values.first-left, len);   
    std::string substcol = col.substr(values.second-left, len); 
    
    if(len >= minLEN)
    {
        char *read1 = new char[len+1];
        char *read2 = new char[len+1];
        strcpy(read1, substrow.c_str());
        strcpy(read2, substcol.c_str());
    
        /* In len we expect len*pow(pCORR, 2) correct base-pairs, so here we compute the maximum 
        edit distance as len - mean + 2 standard deviation = len - (len*pow(pCORR, 2) + 2 * qrt(len*pow(pCORR, 2)*(1-pow(pCORR, 2)) */
        int maxed = len-len*pCORR+2*sqrt(len*pCORR*pWRON); 
    
        EdlibAlignResult result = edlibAlign(read1, len, read2, len, edlibNewAlignConfig(maxed, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
        delete [] read1;
        delete [] read2;
    
        if(result.editDistance != -1)
        {
            align = true;
            ed = result.editDistance;
        }
        else align = false;
    
        edlibFreeAlignResult(result);
    }
    else align = false;

    return align;
}
#endif

/* CSV containing CSC indices of the output sparse matrix*/
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
void HeapSpGEMM(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, std::vector<string> & reads, FT & getvaluetype)
{   
    #ifdef _OSX /* OSX-based memory consumption implementation */
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
    
    #ifdef _LINUX /* LINUX-based memory consumption implementation */
    if(sysinfo(&info) != 0)
    {
        return false;
    }   
    unsigned long free_memory = info.freeram * info.mem_unit;
    free_memory += info.freeswap * info.mem_unit;
    free_memory += info.bufferram * info.mem_unit;
    #endif

    uint64_t d = B.nnz/B.cols; // about d nnz each col
    uint64_t rsv = B.nnz*d;    // worst case
    uint64_t required_memory = (rsv)*sizeof(size_t);
    IT * rowids;
    FT * values;
    
    if(required_memory > free_memory)
    {
        
        int numBlocks = required_memory/free_memory;
        int colsPerBlock = B.cols/numBlocks;                // define number of columns for each blocks

        /* multi thread variable definition */
        int * colStart = new int[numBlocks+1];              // need one block more for remaining cols
        int * colEnd = new int[numBlocks+1];                // need one block more for remaining cols
        int * numCols = new int[numBlocks+1];

        cout << "numBlocks " << numBlocks+1 << endl;
        cout << "colsPerBlock " << colsPerBlock << endl;

        colStart[0] = 0;
        colEnd[0] = 0;
        numCols[0] = 0;

        int colsTrace = 0;
        for(int i = 0; i < numBlocks+1; ++i)
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

        //IT * colptr = new IT[B.cols+1]; 
        //colptr[0] = 0;
        std::stringstream myBatch;

        #ifdef _MULTPTR
        std::stringstream myTrue;
        std::stringstream myFalse;
        #endif

        #pragma omp parallel for private(myBatch) shared(colStart,colEnd,numCols)
        for(int b = 0; b < numBlocks+1; ++b) 
        { 
            vector<IT> * RowIdsofC = new vector<IT>[numCols[b]];  // row ids for each column of C (bunch of cols)
            vector<FT> * ValuesofC = new vector<FT>[numCols[b]];  // values for each column of C (bunch of cols)
            
            IT * colptr = new IT[numCols[b]+1];
            colptr[0] = 0;
           
            LocalSpGEMM(colStart[b], colEnd[b], numCols[b], A, B, multop, addop, RowIdsofC, ValuesofC);
     
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

            delete [] RowIdsofC;
            delete [] ValuesofC;

            /* Local Alignment before write on stringstream */
            for(int i=0; i<numCols[b]; ++i) {
                for(int j=colptr[i]; j<colptr[i+1]; ++j) {
                    #ifdef _EDLIB
                    #ifdef _MULTPTR
                    int len, ed;
                    if(edlibOp(values[j]->front().second, reads[rowids[j]], reads[i+colStart[b]], len, ed) == true)
                    {
                        if(values[j]->size() > 1)
                        {
                        double ratio = getRatio(values[j]);
                        myTrue << i+colStart[b] << ',' << rowids[j] << ',' << ratio << endl;
                        }
                        myBatch << i+colStart[b] << ',' << rowids[j] << ',' << values[j]->size() << endl;

                    } else if(edlibOp(values[j]->front().second, reads[rowids[j]], reads[i+colStart[b]], len, ed) == false)
                    {
                        if(values[j]->size() > 1)
                        {
                            double ratio = getRatio(values[j]);
                            myFalse << i+colStart[b] << ',' << rowids[j] << ',' << ratio << endl;
                        }
                    }
                    #else
                    if(values[j].first < 2)
                    {
                        int len, ed; 
                        if(edlibOp(values[j].second, reads[rowids[j]], reads[i+colStart[b]], len, ed) == true)
                            myBatch << i+colStart[b] << ',' << rowids[j] << ',' << values[j].first << endl;
                    } 
                    else myBatch << i+colStart[b] << ',' << rowids[j] << ',' << values[j].first << endl;
                    #endif
                    #else
                    // add SeqAn seed and extend alignment
                    if(values[j]->count < 2)
                    {
                        seqan::CharString col = reads[i+colStart[b]];
                        seqan::CharString row = reads[rowids[j]];
                        // Seed creation (as of now, we use just one seed)
                        seqan::Seed<seqan::Simple> seed(values[j]->pos[0], values[j]->pos[1], KMER_LENGTH); // begin(read i), begin(read j) and length(seed)
                        // void extendSeed(seed, database, query, direction, MatchExtend);
                        // void extendSeed(seed, database, query, direction, scoreMatrix, scoreDropOff, {UngappedXDrop, GappedXDrop});
                        seqan::Score<int, seqan::Simple> scoringScheme(1, -1, -1);
                        // Select drop off: max(\epsilon, k * error_rate * sequence_extended_so_far)
                        // Some initial guesses would be \epsilon=3, and k=2 (these are of course tunable)
                        extendSeed(seed, col, row, seqan::EXTEND_BOTH, scoringScheme, 90, seqan::UnGappedXDrop());
                        if(max(endPositionH(seed)-beginPositionH(seed), endPositionV(seed)-beginPositionV(seed)) > 300)
                            myBatch << i+colStart[b] << ',' << rowids[j] << ',' << values[j]->count << ',' << values[j]->pos[0] << ',' << values[j]->pos[1] << ',' << values[j]->pos[2] << ',' << values[j]->pos[3] << endl;
                        // std::cout << endPositionH(seed3) << std::endl;  //output: 14
                        // std::cout << endPositionV(seed3) << std::endl;  //output: 13
                        // Comment from SeqAn developer: for ungapped X-drop extension of Simple Seeds, we can simply
                        // update the begin and end values in each dimension. TODO: understand how to discriminare "valide" reads pairs.
                        // TODO: re-define macros for bits
                       
                    }
                    else 
                    myBatch << i+colStart[b] << ',' << rowids[j] << ',' << values[j]->count << ',' << values[j]->pos[0] << ',' << values[j]->pos[1] << ',' << values[j]->pos[2] << ',' << values[j]->pos[3] << endl;
                    #endif
                }
            }

            delete [] colptr;
            delete [] rowids;
            delete [] values;

            #pragma omp critical
            {
                #ifdef _MULTPTR
                writeToFile(myFalse, "falsesample.csv");
                writeToFile(myTrue, "truesample.csv");
                myFalse.str(std::string());
                myTrue.str(std::string());
                #endif
                writeToFile(myBatch, "spmat.csv");
                myBatch.str(std::string());
            }
        }
        //delete [] colptr;
        delete [] colStart;
        delete [] colEnd;
        delete [] numCols;
    } else
    {
        int colStart = 0;
        int colEnd = B.cols;
        int numCols = B.cols;
        vector<IT> * RowIdsofC = new vector<IT>[B.cols];      // row ids for each column of C
        vector<FT> * ValuesofC = new vector<FT>[B.cols];      // values for each column of C

        LocalSpGEMM(colStart, colEnd, numCols, A, B, multop, addop, RowIdsofC, ValuesofC);

        IT * colptr = new IT[B.cols+1];
        colptr[0] = 0;
        std::stringstream myBatch;

        for(int i=0; i < B.cols; ++i)        // for all edge lists (do NOT parallelize)
            colptr[i+1] = colptr[i] + RowIdsofC[i].size();

        IT * rowids = new IT[colptr[B.cols]];
        FT * values = new FT[colptr[B.cols]];
        
        for(int i=0; i < B.cols; ++i)        // for all edge lists (do NOT parallelize)
        {
            copy(RowIdsofC[i].begin(), RowIdsofC[i].end(), rowids + colptr[i]);
            copy(ValuesofC[i].begin(), ValuesofC[i].end(), values + colptr[i]);
        }

        delete [] RowIdsofC;
        delete [] ValuesofC;

        for(int i=0; i<B.cols; ++i) {
            for(int j=colptr[i]; j<colptr[i+1]; ++j) {
                #ifdef _EDLIB
                #ifdef _MULTPTR
                if(values[j]->size() < 2)
                {
                    int len, ed;
                    if(edlibOp(values[j]->front().second, reads[rowids[j]], reads[i], len, ed) == true)
                    myBatch << i << ',' << rowids[j] << ',' << values[j]->size() << endl;
                }
                else myBatch << i << ',' << rowids[j] << ',' << values[j]->size() << endl;
                #else
                if(values[j].first < 2)
                {
                int len, ed; 
                if(edlibOp(values[j].second, reads[rowids[j]], reads[i], len, ed) == true)
                    myBatch << i << ',' << rowids[j] << ',' << values[j].first << endl;
                }
                else myBatch << i << ',' << rowids[j] << ',' << values[j].first << endl;
                #endif
                #else
                myBatch << i << ',' << rowids[j] << ',' << values[j]->count << ',' << values[j]->pos[0] << ',' << values[j]->pos[1] << ',' << values[j]->pos[2] << ',' << values[j]->pos[3] << endl;
                #endif
            }
        }
        
        delete [] rowids;
        delete [] values;
        delete [] colptr;

        writeToFile(myBatch, "spmat.csv");
        myBatch.str(std::string());
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