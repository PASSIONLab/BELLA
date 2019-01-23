#include "CSC.h"
#include "align.h"
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

#ifndef PRINT
#define PRINT
#endif

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

double safety_net = 1.2;

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
	int ithread = omp_get_thread_num();

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
IT* estimateFLOP(const CSC<IT,NT> & A, const CSC<IT,NT> & B, bool lowtriout)
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
    
	IT* colflopC = new IT[B.cols]; // nnz in every  column of C
    	
	#pragma omp parallel for
   	for(IT i=0; i< B.cols; ++i)
    	{
        	colflopC[i] = 0;
    	}

	#pragma omp parallel for
    	for(IT i=0; i < B.cols; ++i)
    	{
        	size_t nnzcolB = B.colptr[i+1] - B.colptr[i]; //nnz in the current column of B
		int myThread = omp_get_thread_num();
		for (IT j = B.colptr[i]; j < B.colptr[i+1]; ++j)	// all nonzeros in that column of B
		{
			IT col2fetch = B.rowids[j];	// find the row index of that nonzero in B, which is the column to fetch in A
			IT nnzcolA = 0;

			if(lowtriout)
			{
				for(IT k = A.colptr[col2fetch]; k < A.colptr[col2fetch+1]; ++k) // all nonzeros in this column of A
				{
					// i is the column_id of the output and A.rowids[k] is the row_id of the output
					if(i < A.rowids[k])
					{
						++nnzcolA;
					}
				}
			}
			else
			{
				nnzcolA =  A.colptr[col2fetch+1]- A.colptr[col2fetch]; // nonzero count of that column of A
			}
			colflopC[i] += nnzcolA;
		}
	}
    	return colflopC;
}

// estimate space for result of SpGEMM with Hash
template <typename IT, typename NT>
IT* estimateNNZ_Hash(const CSC<IT,NT> & A, const CSC<IT,NT> & B, const size_t *flopC, bool lowtriout)
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

        for(size_t j=0; j < ht_size; ++j)
        {
            globalHashVec[j] = -1;
        }
            
	for (IT j = B.colptr[i]; j < B.colptr[i+1]; ++j)	// all nonzeros in that column of B
	{
		IT col2fetch = B.rowids[j];	// find the row index of that nonzero in B, which is the column to fetch in A
		for(IT k = A.colptr[col2fetch]; k < A.colptr[col2fetch+1]; ++k) // all nonzeros in this column of A
		{
			IT key = A.rowids[k];

			if(lowtriout && i >= key)	// i is the column_id of the output and key is the row_id of the output
				continue;

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

//! Hash based column-by-column spgemm algorithm. Based on earlier code by Buluc, Azad, and Nagasaka
//! If lowtriout= true, then only creates the lower triangular part: no diagonal and no upper triangular
template <typename IT, typename NT, typename MultiplyOperation, typename AddOperation, typename FT>
void LocalSpGEMM(IT & start, IT & end, const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, 
		vector<IT> * RowIdsofC, vector<FT> * ValuesofC, IT* colptrC, bool lowtriout)
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
	std::vector< std::pair<IT,FT>> globalHashVec(ht_size);
           
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

			if(lowtriout && i >= key)	// i is the column_id of the output and key is the row_id of the output
				continue;

			FT result =  multop(A.values[k], valueofB);
                	IT hash = (key*hashScale) & (ht_size-1);
                	while (1) //hash probing
                	{
                    		if (globalHashVec[hash].first == key) //key is found in hash table
                    		{
                            		globalHashVec[hash].second = addop(result, globalHashVec[hash].second);
                        		break;
                    		}
                    		else if (globalHashVec[hash].first == -1) //key is not registered yet
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
       // gather non-zero elements from hash table (and then sort them by row indices if needed)
       IT index = 0;
       for (IT j=0; j < ht_size; ++j)
       {
	       if (globalHashVec[j].first != -1)
	       {
		       globalHashVec[index++] = globalHashVec[j];
	       }
       }
#ifdef SORTCOLS
       std::sort(globalHashVec.begin(), globalHashVec.begin() + index, sort_less<IT, NT>);
#endif
       RowIdsofC[i-start].resize(index); 
       ValuesofC[i-start].resize(index);
 
       for (IT j=0; j< index; ++j)
       {
	       RowIdsofC[i-start][j] = globalHashVec[j].first;
	       ValuesofC[i-start][j] = globalHashVec[j].second;
       }
    }

}

double estimateMemory(const BELLApars & b_pars)
{
    double free_memory;
    if (b_pars.userDefMem)
    {
    	free_memory = b_pars.totalMemory * 1024 * 1024;
    }
    else
    {
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
        free_memory = (double) vm_stats.free_count * (double)page_size;
    } 
#elif defined (LINUX) // LINUX-based memory consumption implementation 
    if(sysinfo(&info) != 0)
    {
        return false;
    }   
    free_memory = info.freeram * info.mem_unit;
    free_memory += info.freeswap * info.mem_unit;
    free_memory += info.bufferram * info.mem_unit;
#else
    free_memory = b_pars.totalMemory * 1024 * 1024;	// memory is neither user-supplied nor can be estimated, so use BELLA's default
#endif
    }
    return free_memory;
}



void PostAlignDecision(const seqAnResult & maxExtScore, const readType_ & read1, const readType_ & read2, 
					const BELLApars & b_pars, double ratioPhi, int count, stringstream & myBatch, size_t & outputted,
					size_t & numBasesAlignedTrue, size_t & numBasesAlignedFalse)
{
	auto maxseed = maxExtScore.seed;	// returns a seqan:Seed object

	// {begin/end}Position{V/H}: Returns the begin/end position of the seed in the query (vertical/horizonral direction)
	// these four return seqan:Tposition objects
	auto begpV = beginPositionV(maxseed);
	auto endpV = endPositionV(maxseed);	
	auto begpH = beginPositionH(maxseed);
	auto endpH = endPositionH(maxseed);

	// get references for better naming
	const string& seq1 = read1.seq;
	const string& seq2 = read2.seq;
			
	int read1len = seq1.length();
	int read2len = seq2.length();	
	bool passed = false;

	if(b_pars.adapThr)
	{
		int diffCol = endpV - begpV;
		int diffRow = endpH - begpH;
		int minLeft = min(begpV, begpH);
		int minRight = min(read2len - endpV, read1len - endpH);

		int ov = minLeft+minRight+(diffCol+diffRow)/2;
		double newThr = (1-b_pars.deltaChernoff)*(ratioPhi*(double)ov);

		if((double)maxExtScore.score > newThr)
		{
			if(b_pars.alignEnd)
			{
				if(toEnd(begpV, endpV, read2len, begpH, endpH, read1len, b_pars.relaxMargin))
					passed = true;		
			}
			else 
			{
				passed = true;
			}
		}
	}
	else if(maxExtScore.score > b_pars.defaultThr)
	{
		if(b_pars.alignEnd)
		{
			if(toEnd(begpV, endpV, read2len, begpH, endpH, read1len, b_pars.relaxMargin))
				passed = true;
		}
		else 
		{
			passed = true;
		}
	}
	if(passed)
	{
		myBatch << read2.nametag << '\t' << read1.nametag << '\t' << count << '\t' << maxExtScore.score << '\t' << maxExtScore.strand << '\t' << 
			begpV << '\t' << endpV << '\t' << read2len << '\t' << begpH << '\t' << endpH << '\t' << read1len << endl;
		++outputted;
		numBasesAlignedTrue += (endpV-begpV);	
	}		
	else
	{
		numBasesAlignedFalse += (endpV-begpV);		
	}		
}

template <typename IT, typename FT>
auto RunPairWiseAlignments(IT start, IT end, IT offset, IT * colptrC, IT * rowids, FT * values, const readVector_ & reads, 
								int kmer_len, int xdrop, char* filename, const BELLApars & b_pars, double ratioPhi)
{
    size_t alignedpairs = 0;
    size_t alignedbases = 0;
    size_t totalreadlen = 0;
    size_t totaloutputt = 0;
    size_t totsuccbases = 0;
    size_t totfailbases = 0;
    
    int numThreads = 1;
#pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }

    vector<stringstream> vss(numThreads);	// any chance of false sharing here? depends on how stringstream is implemented. optimize later if needed...

#pragma omp parallel for	
    for(IT j = start; j<end; ++j) // for (end-start) columns of A^T A (one block)
    {
    	size_t numAlignmentsThread = 0;
        size_t numBasesAlignedThread = 0;
        size_t readLengthsThread = 0;
	size_t numBasesAlignedTrue = 0;
	size_t numBasesAlignedFalse = 0;
	
	size_t outputted = 0;

	int ithread = omp_get_thread_num();	

	for (IT i = colptrC[j]; i < colptrC[j+1]; ++i)	// all nonzeros in that column of A^T A
	{
		size_t rid = rowids[i-offset];	// row id
		size_t cid = j;			// column id
		const string& seq1 = reads[rid].seq;	// get reference for readibility
		const string& seq2 = reads[cid].seq;	// get reference for readibility
			
		int seq1len = seq1.length();
		int seq2len = seq2.length();

		spmatPtr_ val = values[i-offset];		


		if(!b_pars.skipAlignment) // fix -z to not print 
                {	
#ifdef TIMESTEP	
                    	numAlignmentsThread++;
                    	readLengthsThread = readLengthsThread + seq1len + seq2len;
#endif

		    	seqAnResult maxExtScore;
		    	if(val->count == 1)
		    	{
				maxExtScore = seqanAlOne(seq1, seq2, seq1len, val->pos[0], val->pos[1], xdrop, kmer_len);
		    	} 
			else
		    	{
                        	maxExtScore = seqanAlGen(seq1, seq2, seq1len, val->pos[0], val->pos[1], val->pos[2], val->pos[3], xdrop, kmer_len);
		    	}			
#ifdef TIMESTEP
			numBasesAlignedThread += endPositionV(maxExtScore.seed)-beginPositionV(maxExtScore.seed);
#endif
                       	PostAlignDecision(maxExtScore, reads[rid], reads[cid], b_pars, ratioPhi, val->count, vss[ithread], outputted, numBasesAlignedTrue, numBasesAlignedFalse);

		}
	        else 	// if skipAlignment == false do alignment, else save just some info on the pair to file 		
		{
			vss[ithread] << reads[cid].nametag << '\t' << reads[rid].nametag << '\t' << val->count << '\t' << 
                        	seq2len << '\t' << seq1len << endl;
			++outputted;
		}
	} // all nonzeros in that column of A^T A

#ifdef TIMESTEP	
	#pragma omp critical
	{	
        	alignedpairs += numAlignmentsThread;
		alignedbases += numBasesAlignedThread;	
		totalreadlen += readLengthsThread;
		totaloutputt += outputted;
		totsuccbases += numBasesAlignedTrue;		
		totfailbases += numBasesAlignedFalse;
	}
#endif			
    }	// all columns from start...end (omp for loop)


    int64_t * bytes = new int64_t[numThreads];
    for(int i=0; i< numThreads; ++i)
    {
	    vss[i].seekg(0, ios::end);
	    bytes[i] = vss[i].tellg();
	    vss[i].seekg(0, ios::beg);
    }
    int64_t bytestotal = std::accumulate(bytes, bytes+numThreads, static_cast<int64_t>(0));

    std::ofstream ofs(filename, std::ios::binary | std::ios::app);
    cout << "Creating or appending to output file with " << bytestotal << " bytes" << endl;
    ofs.seekp(bytestotal - 1);
    ofs.write("", 1);	// this will likely create a sparse file so the actual disks won't spin yet
    ofs.close();
    
    #pragma omp parallel
    {
	int ithread = omp_get_thread_num();	
	    
    	FILE *ffinal;
	if ((ffinal = fopen(filename, "rb+")) == NULL)	// then everyone fills it
        {
		fprintf(stderr, "File %s failed to open at thread %d\n", filename, ithread);
       	}
	int64_t bytesuntil = std::accumulate(bytes, bytes+ithread, static_cast<int64_t>(0));
	fseek (ffinal , bytesuntil , SEEK_SET );
	std::string text = vss[ithread].str();
	fwrite(text.c_str(),1, bytes[ithread] ,ffinal);
	fflush(ffinal);
	fclose(ffinal);
    }
    delete [] bytes;    
	
    return make_tuple(alignedpairs, alignedbases, totalreadlen, totaloutputt, totsuccbases, totfailbases);
}


/**
  * Sparse multithreaded GEMM.
 **/
template <typename IT, typename NT, typename FT, typename MultiplyOperation, typename AddOperation>
void HashSpGEMM(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, const readVector_ & reads, 
    FT & getvaluetype, int kmer_len, int xdrop, char* filename, const BELLApars & b_pars, double ratioPhi)
{
    double free_memory = estimateMemory(b_pars);

#ifdef PRINT
    cout << "Available RAM is assumed to be : " << free_memory / (1024 * 1024) << " MB" << endl;
#endif

    int numThreads = 1;
#pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }

    IT* flopC = estimateFLOP(A, B, true);
    IT* flopptr = prefixsum<IT>(flopC, B.cols, numThreads);
    IT flops = flopptr[B.cols];

#ifdef PRINT    
    cout << "FLOPS is " << flops << endl;
#endif

    IT* colnnzC = estimateNNZ_Hash(A, B, flopC, true);
    IT* colptrC = prefixsum<IT>(colnnzC, B.cols, numThreads);	// colptrC[i] = rolling sum of nonzeros in C[1...i]
    delete [] colnnzC;
    delete [] flopptr;
    delete [] flopC;
    IT nnzc = colptrC[B.cols];
    double compression_ratio = (double)flops / nnzc;


    uint64_t required_memory = safety_net * nnzc * (sizeof(FT)+sizeof(IT));	// required memory to form the output
    int stages = std::ceil((double) required_memory/ free_memory); 	// form output in stages 
    uint64_t nnzcperstage = free_memory / (safety_net * (sizeof(FT)+sizeof(IT)));

#ifdef PRINT
    cout << "nnz(output): " << nnzc << ", free_memory: " << free_memory << ", required_memory: " << required_memory << endl; 
    cout << "stages: " << stages << ", max nnz per stage: " << nnzcperstage << endl;    
#endif


    IT * colStart = new IT[stages+1];	// one array is enough to set stage boundaries	              
    colStart[0] = 0;

    for(int i = 1; i < stages; ++i)	// colsPerStage is no longer fixed (helps with potential load imbalance)
    {
	// std::upper_bound returns an iterator pointing to the first element 
	// in the range [first, last) that is greater than value, or last if no such element is found
	auto upper = std::upper_bound(colptrC, colptrC+B.cols+1, i*nnzcperstage ); 
	colStart[i]  = upper - colptrC - 1;	// we don't want the element that exceeds our budget, we want the one just before that
    }
    colStart[stages] = B.cols;


    for(int b = 0; b < stages; ++b) 
    {
#ifdef TIMESTEP
        double ovl = omp_get_wtime();
#endif
        vector<IT> * RowIdsofC = new vector<IT>[colStart[b+1]-colStart[b]];    // row ids for each column of C (bunch of cols)
        vector<FT> * ValuesofC = new vector<FT>[colStart[b+1]-colStart[b]];    // values for each column of C (bunch of cols)

        LocalSpGEMM(colStart[b], colStart[b+1], A, B, multop, addop, RowIdsofC, ValuesofC, colptrC, true);

#ifdef TIMESTEP
        double ov2 = omp_get_wtime();
	cout << "Columns [" << colStart[b] << " - " << colStart[b+1] << "] overlap time: " << ov2-ovl << "s" << endl;
#endif
        
	IT endnz = colptrC[colStart[b+1]];
	IT begnz = colptrC[colStart[b]];

        IT * rowids = new IT[endnz-begnz];
        FT * values = new FT[endnz-begnz];

        for(IT i=colStart[b]; i<colStart[b+1]; ++i) // combine step
        {
	    IT loccol = i-colStart[b];
	    IT locnz = colptrC[i]-begnz;
            copy(RowIdsofC[loccol].begin(), RowIdsofC[loccol].end(), rowids + locnz);
            copy(ValuesofC[loccol].begin(), ValuesofC[loccol].end(), values + locnz);
	}

        delete [] RowIdsofC;
        delete [] ValuesofC;

	tuple<size_t, size_t, size_t, size_t, size_t, size_t> alignstats; // (alignedpairs, alignedbases, totalreadlen, outputted, alignedtrue, alignedfalse)
	alignstats = RunPairWiseAlignments(colStart[b], colStart[b+1], begnz, colptrC, rowids, values, reads, kmer_len, xdrop, filename, b_pars, ratioPhi);

#ifdef TIMESTEP
        if(!b_pars.skipAlignment)
        {
	    	double elapsed = omp_get_wtime()-ov2;
            	cout << "Columns [" << colStart[b] << " - " << colStart[b+1] << "] alignment time: " << elapsed << " s | alignment rate: " << static_cast<double>(get<1>(alignstats))/elapsed;
	   	cout << " bases/s | average read length: " <<static_cast<double>(get<2>(alignstats))/(2* get<0>(alignstats));
		cout << " | read pairs aligned this stage: " << get<0>(alignstats) << endl;
		cout << "Outputted " << get<3>(alignstats) << " lines" << endl;
		cout << "Average length of successful alignment " << static_cast<double>(get<4>(alignstats)) / get<3>(alignstats) << " bps" << endl;
		cout << "Average length of failed alignment " << static_cast<double>(get<5>(alignstats)) / (get<0>(alignstats) - get<3>(alignstats)) << " bps" << endl;		
	}
#endif

        delete [] rowids;
        delete [] values;

    }//for(int b = 0; b < states; ++b)

    delete [] colptrC;
    delete [] colStart;
}
