#include "CSC.h"
#include "align.h"
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
#define TIMESTEP
//#define PRINT
//#define THREADLIMIT
//#define MAX_NUM_THREAD 1
//#define OSX
//#define LINUX
//#define RAM

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

extern int totalMemory;


/*
 Multithreaded prefix sum
 Inputs:
    in: an input array
    size: the length of the input array "in"
    nthreads: number of threads used to compute the prefix sum
 
 Output:
    return an array of size "size+1"
    the memory of the output array is allocated internallay
 
 Example:
 
    in = [2, 1, 3, 5]
    out = [0, 2, 3, 6, 11]
 */
template <typename T>
T* prefixsum(T* in, int size, int nthreads)
{
    std::vector<T> tsum(nthreads+1);
    tsum[0] = 0;
    T* out = new T[size+1];
    out[0] = 0;
    T* psum = &out[1];
#pragma omp parallel
    {
	int omp_get_thread_num();

        T sum = 0;
#pragma omp for schedule(static)
        for (int i=0; i<size; i++)
        {
            sum += in[i];
            psum[i] = sum;
        }
        
        tsum[ithread+1] = sum;
#pragma omp barrier
        T offset = 0;
        for(int i=0; i<(ithread+1); i++)
        {
            offset += tsum[i];
        }
		
#pragma omp for schedule(static)
        for (int i=0; i<size; i++)
        {
            psum[i] += offset;
        }
    
    }
    return out;
}

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

// estimate the number of floating point operations of SpGEMM
template <typename IT, typename NT>
IT* estimateFLOP(const CSC<IT,NT> & A, const CSC<IT,NT> & B)
{
	if(A.isEmpty() || B.isEmpty())
	{
		return NULL;
    	}
    	
	int numThreads = 1;
	#ifdef THREADED
	#pragma omp parallel
    	{
        	numThreads = omp_get_num_threads();
   	}
	#endif
    
	IT* colflopC = new IT[B.cols]; // nnz in every  column of C
    	
	#ifdef THREADED
	#pragma omp parallel for
	#endif
   	for(IT i=0; i< B.cols; ++i)
    	{
        	colflopC[i] = 0;
    	}

	#ifdef THREADED
	#pragma omp parallel for
	#endif
    	for(IT i=0; i < B.cols; ++i)
    	{
        	size_t nnzcolB = B.colptr[i+1] - B.colptr[i]; //nnz in the current column of B
		int myThread = 0;
#ifdef THREADED
        	myThread = omp_get_thread_num();
#endif
		for (IT j = B.colptr[i]; j < B.colptr[i+1]; ++j)	// all nonzeros in that column of B
		{
			IT col2fetch = B.rowids[j];	// find the row index of that nonzero in B, which is the column to fetch in A
			IT nnzcolA =  A.colptr[col2fetch+1]- A.colptr[col2fetch]; // nonzero count of that column of A
			colflopC[i] += nnzcolA;
		}
	}
    }
    return colflopC;
}

// estimate space for result of SpGEMM with Hash
template <typename IT, typename NT>
IT* estimateNNZ_Hash(const CSC<IT,NT> & A, const CSC<IT,NT> & B, const size_t *flopC)
{
    if(A.isEmpty() || B.isEmpty())
    {
        return NULL;
    }
    		
    int numThreads = 1;
#pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }

    IT* colnnzC = new IT[B.cols]; // nnz in every  column of C
	
#pragma omp parallel for
    for(IT i=0; i< B.cols; ++i)
    {
        colnnzC[i] = 0;
    } 

#pragma omp parallel for
    for(IT i=0; i < B.cols; ++i)	// for each column of B
    {
        size_t nnzcolB = B.colptr[i+1] - B.colptr[i]; //nnz in the current column of B
	int myThread = omp_get_thread_num();
    		
        // Hash
        const size_t minHashTableSize = 16;
        const size_t hashScale = 107;

        // Initialize hash tables
        size_t ht_size = minHashTableSize;
        while(ht_size < flopC[i]) //ht_size is set as 2^n
        {
            ht_size <<= 1;
        }
        std::vector<IT> globalHashVec(ht_size);

        for(size j=0; j < ht_size; ++j)
        {
            globalHashVec[j] = -1;
        }
            
	for (IT j = B.colptr[i]; j < B.colptr[i+1]; ++j)	// all nonzeros in that column of B
	{
		IT col2fetch = B.rowids[j];	// find the row index of that nonzero in B, which is the column to fetch in A
		for(IT k = A.colptr[col2fetch]; k < A.colptr[col2fetch+1]; ++k) // all nonzeros in this column of A
		{
			IT key = A.rowids[k];
                	IT hash = (key*hashScale) & (ht_size-1);
                	while (1) //hash probing
                	{
                    		if (globalHashVec[hash] == key) //key is found in hash table
                    		{
                        		break;
                    		}
                    		else if (globalHashVec[hash] == -1) //key is not registered yet
                    		{
                        		globalHashVec[hash] = key;
                        		colnnzC[i] ++;
                        		break;
                    		}
                    		else //key is not found
                    		{
                        		hash = (hash+1) & (ht_size-1);	// don't exit the while loop yet
                    		}
			}
		}
	}
    }
    
    return colnnzC;
}

template <typename IT, typename NT, typename MultiplyOperation, typename AddOperation, typename FT>
void LocalSpGEMM(IT & start, IT & end, const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, 
		vector<IT> * RowIdsofC, vector<FT> * ValuesofC, IT* colptrC)
{
#pragma omp parallel for	
    for(IT i = start; i<end; ++i) // for bcols of B (one block)
    {
        const IT minHashTableSize = 16;
	const IT hashScale = 107;
    	size_t nnzcolC = colptrC[i+1] - colptrC[i]; //nnz in the current column of C (=Output)
            
	IT ht_size = minHashTableSize;
	while(ht_size < nnzcolC) //ht_size is set as 2^n
	{
		ht_size <<= 1;
	}
	std::vector< std::pair<IT,NT>> globalHashVec(ht_size);
           
	// Initialize hash tables
	for(IT j=0; j < ht_size; ++j)
	{
		globalHashVec[j].first = -1;
	}
            
	for (IT j = B.colptr[i]; j < B.colptr[i+1]; ++j)	// all nonzeros in that column of B
	{
		IT col2fetch = B.rowids[j];	// find the row index of that nonzero in B, which is the column to fetch in A
		NT valueofB = B.values[j];
		for(IT k = A.colptr[col2fetch]; k < A.colptr[col2fetch+1]; ++k) // all nonzeros in this column of A
		{
			IT key = A.rowids[k];
			NT result =  multop(A.values[k], valueofB);
                	IT hash = (key*hashScale) & (ht_size-1);
                	while (1) //hash probing
                	{
                    		if (globalHashVec[hash].first == key) //key is found in hash table
                    		{
                            		globalHashVec[hash].second = addop(result, globalHashVec[hash].second);
                        		break;
                    		}
                    		else if (globalHashVec[hash] == -1) //key is not registered yet
                    		{
                        		globalHashVec[hash].first = key;
                        		globalHashVec[hash].second = result;					
                        		break;
                    		}
                    		else //key is not found
                    		{
                        		hash = (hash+1) & (ht_size-1);	// don't exit the while loop yet
                    		}
			}
		}
	}
       // gather non-zero elements from hash table, and then sort them by row indices
       IT index = 0;
       for (IT j=0; j < ht_size; ++j)
       {
	       if (globalHashVec[j].first != -1)
	       {
		       globalHashVec[index++] = globalHashVec[j];
	       }
       }
       std::sort(globalHashVec.begin(), globalHashVec.begin() + index, sort_less<IT, NT>);
       RowIdsofC[i-start].resize(index); 
       ValuesofC[i-start].resize(index);
 
       for (IT j=0; j< index; ++j)
       {
	       RowIdsofC[i] = globalHashVec[j].first;
	       ValuesofC[i] = globalHashVec[j].second;
       }
    }

}

uint64_t estimateMemory()
{
    uint64_t free_memory;
#if defined (OSX) // OSX-based memory consumption implementation 
    vm_size_t page_size;
    mach_port_t mach_port;
    mach_msg_type_number_t count;
    vm_statistics64_data_t vm_stats;

    mach_port = mach_host_self();
    count = sizeof(vm_stats) / sizeof(natural_t);

    if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
                KERN_SUCCESS == host_statistics64(mach_port, HOST_VM_INFO,
                                                (host_info64_t)&vm_stats, &count))
    {
        free_memory = (int64_t)vm_stats.free_count * (int64_t)page_size;
    } 
#elif defined (LINUX) // LINUX-based memory consumption implementation 
    if(sysinfo(&info) != 0)
    {
        return false;
    }   
    unsigned long free_memory = info.freeram * info.mem_unit;
    free_memory += info.freeswap * info.mem_unit;
    free_memory += info.bufferram * info.mem_unit;
#else
    free_memory = static_cast<uint64_t>(total_memory) * 1024 * 1024;
#endif
    return free_memory;
}


void RunPairWiseAlignments()
{
	for(int i=0; i<numCols[b]; ++i)
        { 
            for(int j=colptr[i]; j<colptr[i+1]; ++j) 
            {
                if(!skipAlignment) // fix -z to not print 

                {
                    /* Local alignment before write on stringstream */
#ifdef TIMESTEP
                    /* Progress report */
                    numAlignmentsThread++;
                    readLengthsThread += globalInstance->at(rowids[j]).seq.length();
                    readLengthsThread += globalInstance->at(i+colStart[b]).seq.length();
#endif
                    if(values[j]->count == 1)
                    {
                        longestExtensionScore = seqanAlOne(globalInstance->at(rowids[j]).seq, globalInstance->at(i+colStart[b]).seq, 
                            globalInstance->at(rowids[j]).seq.length(), values[j]->pos[0], values[j]->pos[1], xdrop, kmer_len);
#ifdef TIMESTEP
                        /* Progress report */
                        numBasesAlignedThread += endPositionV(longestExtensionScore.seed)-beginPositionV(longestExtensionScore.seed);
#endif
                        if(adapThr)
                        {   // "function" this
                            int diffCol = endPositionV(longestExtensionScore.seed)-beginPositionV(longestExtensionScore.seed);
                            int diffRow = endPositionH(longestExtensionScore.seed)-beginPositionH(longestExtensionScore.seed);
                            int minLeft = min(beginPositionV(longestExtensionScore.seed), beginPositionH(longestExtensionScore.seed));
                            int minRight = min(globalInstance->at(i+colStart[b]).seq.length()-endPositionV(longestExtensionScore.seed), globalInstance->at(rowids[j]).seq.length()-endPositionH(longestExtensionScore.seed));

                            int ov = minLeft+minRight+(diffCol+diffRow)/2;
                            double newThr = (1-deltaChernoff)*(ratioPhi*(double)ov); 

                            if((double)longestExtensionScore.score > newThr)
                            {
                                if(alignEnd)
                                {
                                    bool aligntoEnd = toEnd(beginPositionV(longestExtensionScore.seed), endPositionV(longestExtensionScore.seed), globalInstance->at(i+colStart[b]).seq.length(), 
                                        beginPositionH(longestExtensionScore.seed), endPositionH(longestExtensionScore.seed), globalInstance->at(rowids[j]).seq.length(),relaxMargin);
                                
                                    if(aligntoEnd)
                                        myBatch << globalInstance->at(i+colStart[b]).nametag << '\t' << globalInstance->at(rowids[j]).nametag << '\t' << values[j]->count << '\t' << longestExtensionScore.score << '\t' << longestExtensionScore.strand << '\t' << beginPositionV(longestExtensionScore.seed) << '\t' << 
                                            endPositionV(longestExtensionScore.seed) << '\t' << globalInstance->at(i+colStart[b]).seq.length() << '\t' << beginPositionH(longestExtensionScore.seed) << '\t' << endPositionH(longestExtensionScore.seed) <<
                                                '\t' << globalInstance->at(rowids[j]).seq.length() << endl;
                                }
                                else myBatch << globalInstance->at(i+colStart[b]).nametag << '\t' << globalInstance->at(rowids[j]).nametag << '\t' << values[j]->count << '\t' << longestExtensionScore.score << '\t' << longestExtensionScore.strand << '\t' << beginPositionV(longestExtensionScore.seed) << '\t' << 
                                        endPositionV(longestExtensionScore.seed) << '\t' << globalInstance->at(i+colStart[b]).seq.length() << '\t' << beginPositionH(longestExtensionScore.seed) << '\t' << endPositionH(longestExtensionScore.seed) <<
                                            '\t' << globalInstance->at(rowids[j]).seq.length() << endl;
                            }
                        }
                        else if(longestExtensionScore.score > defaultThr)
                        {
                            if(alignEnd)
                            {
                                bool aligntoEnd = toEnd(beginPositionV(longestExtensionScore.seed), endPositionV(longestExtensionScore.seed), globalInstance->at(i+colStart[b]).seq.length(), 
                                    beginPositionH(longestExtensionScore.seed), endPositionH(longestExtensionScore.seed), globalInstance->at(rowids[j]).seq.length(),relaxMargin);

                                if(aligntoEnd)
                                    myBatch << globalInstance->at(i+colStart[b]).nametag << '\t' << globalInstance->at(rowids[j]).nametag << '\t' << values[j]->count << '\t' << longestExtensionScore.score << '\t' << longestExtensionScore.strand << '\t' << beginPositionV(longestExtensionScore.seed) << '\t' << 
                                        endPositionV(longestExtensionScore.seed) << '\t' << globalInstance->at(i+colStart[b]).seq.length() << '\t' << beginPositionH(longestExtensionScore.seed) << '\t' << endPositionH(longestExtensionScore.seed) <<
                                            '\t' << globalInstance->at(rowids[j]).seq.length() << endl;
                            }
                            else myBatch << globalInstance->at(i+colStart[b]).nametag << '\t' << globalInstance->at(rowids[j]).nametag << '\t' << values[j]->count << '\t' << longestExtensionScore.score << '\t' << longestExtensionScore.strand << '\t' << beginPositionV(longestExtensionScore.seed) << '\t' << 
                                endPositionV(longestExtensionScore.seed) << '\t' << globalInstance->at(i+colStart[b]).seq.length() << '\t' << beginPositionH(longestExtensionScore.seed) << '\t' << endPositionH(longestExtensionScore.seed) <<
                                    '\t' << globalInstance->at(rowids[j]).seq.length() << endl;
                        }
                    } // if(values[j]->count == 1)
                    else
                    {
                        longestExtensionScore = seqanAlGen(globalInstance->at(rowids[j]).seq, globalInstance->at(i+colStart[b]).seq, 
                            globalInstance->at(rowids[j]).seq.length(), values[j]->pos[0], values[j]->pos[1], values[j]->pos[2], values[j]->pos[3], xdrop, kmer_len);

#ifdef TIMESTEP
                        /* Progress report */
                        numBasesAlignedThread += endPositionV(longestExtensionScore.seed)-beginPositionV(longestExtensionScore.seed);
#endif
                        if(adapThr)
                        {    // "function" this
                            int diffCol = endPositionV(longestExtensionScore.seed)-beginPositionV(longestExtensionScore.seed);
                            int diffRow = endPositionH(longestExtensionScore.seed)-beginPositionH(longestExtensionScore.seed);
                            int minLeft = min(beginPositionV(longestExtensionScore.seed), beginPositionH(longestExtensionScore.seed));
                            int minRight = min(globalInstance->at(i+colStart[b]).seq.length()-endPositionV(longestExtensionScore.seed), globalInstance->at(rowids[j]).seq.length()-endPositionH(longestExtensionScore.seed));

                            int ov = minLeft+minRight+(diffCol+diffRow)/2;
                            double newThr = (1-deltaChernoff)*(ratioPhi*(double)ov); 

                            if((double)longestExtensionScore.score > newThr)
                            {
                                myBatch << globalInstance->at(i+colStart[b]).nametag << '\t' << globalInstance->at(rowids[j]).nametag << '\t' << values[j]->count << '\t' << longestExtensionScore.score << '\t' << longestExtensionScore.strand << '\t' << beginPositionV(longestExtensionScore.seed) << '\t' << 
                                        endPositionV(longestExtensionScore.seed) << '\t' << globalInstance->at(i+colStart[b]).seq.length() << '\t' << beginPositionH(longestExtensionScore.seed) << '\t' << endPositionH(longestExtensionScore.seed) <<
                                            '\t' << globalInstance->at(rowids[j]).seq.length() << endl;
                            }
                        }
                        else if(longestExtensionScore.score > defaultThr)
                        {
                            myBatch << globalInstance->at(i+colStart[b]).nametag << '\t' << globalInstance->at(rowids[j]).nametag << '\t' << values[j]->count << '\t' << longestExtensionScore.score << '\t' << longestExtensionScore.strand << '\t' << beginPositionV(longestExtensionScore.seed) << '\t' << 
                                endPositionV(longestExtensionScore.seed) << '\t' << globalInstance->at(i+colStart[b]).seq.length() << '\t' << beginPositionH(longestExtensionScore.seed) << '\t' << endPositionH(longestExtensionScore.seed) <<
                                    '\t' << globalInstance->at(rowids[j]).seq.length() << endl;
                        }
                    }
                }// if skipAlignment == false do alignment, else save just some info on the pair to file
                else myBatch << globalInstance->at(i+colStart[b]).nametag << '\t' << globalInstance->at(rowids[j]).nametag << '\t' << values[j]->count << '\t' << 
                        globalInstance->at(i+colStart[b]).seq.length() << '\t' << globalInstance->at(rowids[j]).seq.length() << endl;
            }// for(int j=colptr[i]; j<colptr[i+1]; ++j)
        }// for(int i=0; i<numCols[b]; ++i)

#ifdef TIMESTEP
        if(!skipAlignment)
        {
#pragma omp critical
        {
            alignSoFar += numAlignmentsThread;
            cout << "[" << omp_get_thread_num()+1 << "] alignment time: " << omp_get_wtime()-ov2 << " s | alignment rate: " << numBasesAlignedThread/(omp_get_wtime()-ov2) << " bases/s | average read length: " << 
                        readLengthsThread/(2*numAlignmentsThread) << " | read pairs aligned so far: " << alignSoFar << endl;
            // rate in pair/s: numAlignmentsThread/(omp_get_wtime()-ov2)
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
}


/**
  * Sparse multithreaded GEMM.
 **/
template <typename IT, typename NT, typename FT, typename MultiplyOperation, typename AddOperation>
void HashSpGEMM(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, readVector_ & read, 
    FT & getvaluetype, int kmer_len, int xdrop, int defaultThr, char* filename, bool skipAlignment, double ratioPhi, bool adapThr, bool alignEnd, int relaxMargin, double deltaChernoff)
{
    size_t upperlimit = 10000000; // in bytes
    int64_t free_memory = EstimateMemory();

#ifdef PRINT
    cout << "Available RAM is assumed to be : " << free_memory / (1024 * 1024) << " MB" << endl;
#endif

    IT* flopC = estimateFLOP(A, B);
    IT* flopptr = prefixsum<IT>(flopC, B.cols, numThreads);
    IT flops = flopptr[B.cols];

#ifdef PRINT    
    cout << "FLOPS is " << flops << endl;
#endif

    IT* colnnzC = estimateNNZ_Hash(A, B, flopC);
    IT* colptrC = prefixsum<IT>(colnnzC, B.cols, numThreads);	// colptrC[i] = rolling sum of nonzeros in C[1...i]
    delete [] colnnzC;
    delete [] flopptr;
    delete [] flopC;
    IT nnzc = colptrC[B.cols];
    double compression_ratio = (double)flops / nnzc;

#ifdef PRINT
    cout << "nnz(output): " << nnzc << endl; 
#endif

    uint64_t required_memory = (nnzc)*sizeof(FT);	// required memory to form the output
    int stages = required_memory/free_memory; 		// form output in stages

    IT * colStart = new IT[stages+1];	// one array is enough to set stage boundaries	              
    colStart[0] = 0;

    for(int i = 1; i < stages; ++i)	// colsPerStage is no longer fixed (helps with potential load imbalance)
    {
	// std::upper_bound returns an iterator to the first element greater than a certain value 
	auto upper = std::upper_bound(colptrC, colptrC+B.cols, i * required_memory); 
	colStart[i]  = upper - colptrC;
    }
    colStart[stages] = B.cols;

    copy(colStart, colStart+stages+1, ostream_iterator<IT>(cout, " "));
    cout << endl;

    shared_ptr<readVector_> globalInstance = make_shared<readVector_>(read); 
    stringstream myBatch; // each thread(?) saves its results in its private stringstream variable
    seqAnResult longestExtensionScore;

    for(int b = 0; b < stages; ++b) 
    {
#ifdef TIMESTEP
        double ovl = omp_get_wtime();
#endif
        vector<IT> * RowIdsofC = new vector<IT>[colStart[b+1]-colStart[b]];    // row ids for each column of C (bunch of cols)
        vector<FT> * ValuesofC = new vector<FT>[colStart[b+1]-colStart[b]];    // values for each column of C (bunch of cols)

        LocalSpGEMM(colStart[b], colStart[b+1], A, B, multop, addop, RowIdsofC, ValuesofC, colptrC);

#ifdef TIMESTEP
        double ov2 = omp_get_wtime();
#pragma omp critical
        {
            cout << "[" << omp_get_thread_num()+1 << "] overlap time: " << ov2-ovl << "s" << endl;
        }
#endif
        
        IT * rowids = new IT[colptr[numCols[b]]];
        FT * values = new FT[colptr[numCols[b]]];

        for(IT i=colStart[b]; i<colStart[b+1]; ++i) // combine step
        {
	    IT locind = i-colStart[b];
            copy(RowIdsofC[locind].begin(), RowIdsofC[locind].end(), rowids + colptrC[i]);
            copy(ValuesofC[locind].begin(), ValuesofC[locind].end(), values + colptrC[i]);
        }

        delete [] RowIdsofC;
        delete [] ValuesofC;

        size_t numAlignmentsThread = 0;
        size_t numBasesAlignedThread = 0;
        size_t readLengthsThread = 0;

	RunPairWiseAlignments();

    }//for(int b = 0; b < states; ++b)

    delete [] colStart;
    delete [] colEnd;
    delete [] numCols; 
}
