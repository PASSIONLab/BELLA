#include "CSC.h"
#include "align.h"
#include "common.h"
#include "../kmercode/hash_funcs.h"
#include "../kmercode/Kmer.hpp"
#include "../kmercode/Buffer.h"
#include "../kmercode/common.h"
#include "../kmercode/fq_reader.h"
#include "../kmercode/ParallelFASTQ.h"
#include "../libcuckoo/cuckoohash_map.hh"
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
#include <set> 

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

double safety_net = 1.5;

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
char complement (char n)
{   
	switch(n)
	{   
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	}   
	assert(false);
	return ' ';
}

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
/* fix according to PAF format */

void toOriginalCoordinates(int& begpH, int& endpH, const int lenH)
{
	unsigned int tmp = begpH;
	begpH = lenH-endpH;
	endpH = lenH-tmp;
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
IT* estimateNNZ_Hash(const CSC<IT,NT>& A, const CSC<IT,NT>& B, const IT* flopC, bool lowtriout)
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
		const unsigned int minHashTableSize = 16;
		const unsigned int hashScale = 107;

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
		for(IT k = A.colptr[col2fetch]; k < A.colptr[col2fetch+1]; ++k)	// all nonzeros in this column of A
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
	for(IT i = start; i < end; ++i)	// for bcols of B (one block)
	{
		const IT minHashTableSize = 16;
		const IT hashScale = 107;
		size_t nnzcolC = colptrC[i+1] - colptrC[i];	//nnz in the current column of C (=Output)

	IT ht_size = minHashTableSize;
	while(ht_size < nnzcolC)	//ht_size is set as 2^n
	{
		ht_size <<= 1;
	}
	std::vector< std::pair<IT,FT>> globalHashVec(ht_size);

	//	Initialize hash tables
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

			//	i is the column_id of the output and key is the row_id of the output
			if(lowtriout && i >= key)
				continue;

			//	GG: modified to get read ids needed to compute alnlenerlap length
			FT result =  multop(A.values[k], valueofB, key, i);

			IT hash = (key*hashScale) & (ht_size-1);
			while (1) //hash probing
			{
				if (globalHashVec[hash].first == key) //key is found in hash table
				{	//	GG: addop temporary modify, remalnlene key, i after testing
					globalHashVec[hash].second = addop(result, globalHashVec[hash].second, key, i);
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
	vm_unsigned int page_size;
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

char
complementbase(char n) {
	switch(n)
	{
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	}
	assert(false);
	return ' ';
}

std::string
reversecomplement(const std::string& seq) {

	std::string cpyseq = seq;
	std::reverse(cpyseq.begin(), cpyseq.end());

	std::transform(
		std::begin(cpyseq),
		std::end  (cpyseq),
		std::begin(cpyseq),
	complementbase);

	return cpyseq;
}


void PostAlignDecisionGPU(const loganResult & maxExtScore, const readType_ & read1, const readType_ & read2, 
					const BELLApars & b_pars, double ratiophi, int count, stringstream & myBatch, size_t & outputted,
					size_t & numBasesAlignedTrue, size_t & numBasesAlignedFalse, bool & passed)
{
	// returns a Logan::Seed object
	SeedL maxseed = maxExtScore.seed;

	// {begin/end}Position{V/H}: Returns the begin/end position of the seed in the query (vertical/horizonral direction)
	// these four return seqan:Tposition objects
	int begpV  	 = getBeginPositionV(maxseed);
	int endpV 	 = getEndPositionV(maxseed);	
	int begpH 	 = getBeginPositionH(maxseed);
	int endpH 	 = getEndPositionH(maxseed);

	// get references for better naming
	const string& seq1 = read1.seq;	// H
	const string& seq2 = read2.seq;	// Vzw

	unsigned short int read1len = seq1.length();
	unsigned short int read2len = seq2.length();

	//	GG: divergence estimation
	unsigned short int overlapLenV = endpV - begpV;
	unsigned short int overlapLenH = endpH - begpH;

	unsigned short int minLeft  = min(begpV, begpH);
	unsigned short int minRight = min(read2len - endpV, read1len - endpH);
	unsigned short int ov       = minLeft + minRight + (overlapLenV + overlapLenH) / 2;

	unsigned short int normLen  = max(overlapLenV, overlapLenH);
	unsigned short int minLen   = min(overlapLenV, overlapLenH);

	if(b_pars.adapThr)
	{
		float mythreshold = (1 - b_pars.deltaChernoff) * (ratiophi * (float)ov);
		if((float)maxExtScore.score > mythreshold)
		{
			passed = true;
		}
	}
	else if(maxExtScore.score > b_pars.defaultThr)
	{
		passed = true;
	}

	if(passed)
	{
		if(!b_pars.outputPaf)		// BELLA output format
		{
			myBatch << read2.nametag << '\t' << read1.nametag << '\t' << count << '\t' << maxExtScore.score << '\t' << ov << '\t' << maxExtScore.strand << '\t' << 
				begpV << '\t' << endpV << '\t' << read2len << '\t' << begpH << '\t' << endpH << '\t' << read1len << endl;
		}
		else
		{
			std::string pafstrand;	// maxExtScore not modifiable
			unsigned short int mapq = 255;			// mapping quality (0-255; 255 for missing)

			if(maxExtScore.strand == "n") pafstrand = "+";
			else pafstrand = "-";

			if(pafstrand == "-")
				toOriginalCoordinates(begpH, endpH, read1len);

			// PAF format is the output format used by minimap/minimap2: https://github.com/lh3/miniasm/blob/master/PAF.md
			myBatch << read2.nametag << '\t' << read2len << '\t' << begpV << '\t' << endpV << '\t' << pafstrand << '\t' << 
				read1.nametag << '\t' << read1len << '\t' << begpH << '\t' << endpH << '\t' << maxExtScore.score << '\t' << ov << '\t' << mapq << endl;
		}
		++outputted;
		numBasesAlignedTrue += (endpV-begpV);
	}
	else
	{
		numBasesAlignedFalse += (endpV-begpV);
	}
}

int getChainLen(const string& seqH, const string& seqV, const int begpH, const int endpH, 
	const int begpV, const int endpV, const int size, const string& strand)
{
	//	GG: TODO modify Kmer to accept different kmerSizes
	cuckoohash_map<Kmer, int> countsubkmers;
	int matchingSubKmers = 0;
	//	std::unordered_mapKmer, bool> countsubkmers;

	for(int i = begpV; i < endpV - size + 1; i++)
	{
		std::string kmerstrfromstr = seqV.substr(i, size);
		Kmer mykmer(kmerstrfromstr.c_str(), kmerstrfromstr.length());
		Kmer lexsmall = mykmer.rep();
		countsubkmers.insert(lexsmall, 1);
	}

	std::string seqHcpy = seqH;
	if(strand == "c")
	{
		seqHcpy = reversecomplement(seqH);
	}

	for(int i = begpH; i < endpH - size + 1; i++)
	{
		std::string kmerstrfromstr = seqHcpy.substr(i, size);
		Kmer mykmer(kmerstrfromstr.c_str(), kmerstrfromstr.length());
		Kmer lexsmall = mykmer.rep();

		int proxy;
		auto found = countsubkmers.find(lexsmall, proxy);
		if(found)
			matchingSubKmers++;
	}

	int chainLen = matchingSubKmers * size;
	return chainLen;
}

void PostAlignDecision(const seqAnResult& maxExtScore, const readType_& read1, const readType_& read2, 
					const BELLApars& b_pars, double ratiophi, int count, stringstream& myBatch, size_t& outputted,
					size_t& numBasesAlignedTrue, size_t& numBasesAlignedFalse, bool& passed, int const& matches)
{
	auto maxseed = maxExtScore.seed;	// returns a seqan:Seed object

	// {begin/end}Position{V/H}: Returns the begin/end position of the seed in the query (vertical/horizonral direction)
	// these four return seqan:Tposition objects
	int begpV = beginPositionV(maxseed);
	int endpV = endPositionV(maxseed);
	int begpH = beginPositionH(maxseed);
	int endpH = endPositionH(maxseed);

	// get references for better naming
	const string& seq1 = read1.seq;	// H
	const string& seq2 = read2.seq;	// Vzw

	unsigned short int read1len = seq1.length();
	unsigned short int read2len = seq2.length();

	//	GG: divergence estimation
	unsigned short int overlapLenV = endpV - begpV;
	unsigned short int overlapLenH = endpH - begpH;

	unsigned short int minLeft  = min(begpV, begpH);
	unsigned short int minRight = min(read2len - endpV, read1len - endpH);
	unsigned short int ov       = minLeft + minRight + (overlapLenV + overlapLenH) / 2;

	unsigned short int normLen  = max(overlapLenV, overlapLenH);
	unsigned short int minLen   = min(overlapLenV, overlapLenH);

	if(b_pars.adapThr)
	{
		float mythreshold = (1 - b_pars.deltaChernoff) * (ratiophi * (float)ov);
		if((float)maxExtScore.score > mythreshold)
		{
			passed = true;
		}
	}
	else if(maxExtScore.score > b_pars.defaultThr)
	{
		passed = true;
	}

	if(passed)
	{
		if(!b_pars.outputPaf)		// BELLA output format
		{
			myBatch << read2.nametag << '\t' << read1.nametag << '\t' << count << '\t' << maxExtScore.score << '\t' << ov << '\t' << maxExtScore.strand << '\t' << 
				begpV << '\t' << endpV << '\t' << read2len << '\t' << begpH << '\t' << endpH << '\t' << read1len << endl;
		}
		else
		{
			std::string pafstrand;	// maxExtScore not modifiable
			unsigned short int mapq = 255;			// mapping quality (0-255; 255 for missing)

			if(maxExtScore.strand == "n") pafstrand = "+";
			else pafstrand = "-";

			if(pafstrand == "-")
				toOriginalCoordinates(begpH, endpH, read1len);

			// PAF format is the output format used by minimap/minimap2: https://github.com/lh3/miniasm/blob/master/PAF.md
			myBatch << read2.nametag << '\t' << read2len << '\t' << begpV << '\t' << endpV << '\t' << pafstrand << '\t' << 
				read1.nametag << '\t' << read1len << '\t' << begpH << '\t' << endpH << '\t' << maxExtScore.score << '\t' << ov << '\t' << mapq << endl;
		}
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
		char* filename, const BELLApars & b_pars, const double& ratiophi)
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

	vector<stringstream> vss(numThreads); // any chance of false sharing here? depends on how stringstream is implemented. optimize later if needed...

#pragma omp parallel for schedule(dynamic)
	for(IT j = start; j < end; ++j)	// for (end-start) columns of A^T A (one block)
	{
		size_t numAlignmentsThread		= 0;
		size_t numBasesAlignedThread	= 0;
		size_t readLengthsThread		= 0;
		size_t numBasesAlignedTrue		= 0;
		size_t numBasesAlignedFalse		= 0;

		size_t outputted = 0;

		int ithread = omp_get_thread_num();

		for (IT i = colptrC[j]; i < colptrC[j+1]; ++i)  // all nonzeros in that column of A^T A
		{
			unsigned int rid = rowids[i-offset];	// row id
			unsigned int cid = j;					// column id

			const string& seq1 = reads[rid].seq;	// get reference for readibility
			const string& seq2 = reads[cid].seq;	// get reference for readibility

			unsigned short int seq1len = seq1.length();
			unsigned short int seq2len = seq2.length();

			spmatPtr_ val = values[i-offset];

			if(!b_pars.skipAlignment) // fix -z to not print 
			{
				numAlignmentsThread++;
				readLengthsThread = readLengthsThread + seq1len + seq2len;

				seqAnResult maxExtScore;
				bool passed = false;

				//	GG: number of matching kmer into the majority voted bin
				unsigned short int matches = val->chain();
				unsigned short int minkmer = 1;

				//	GG: b_pars.minSurvivedKmers should be function of Markov chain
				if (matches < minkmer)
					continue;

				pair<int, int> kmer = val->choose();
				int i = kmer.first, j = kmer.second;

				//	GG: nucleotide alignment
				maxExtScore = alignSeqAn(seq1, seq2, seq1len, i, j, b_pars.xDrop, b_pars.kmerSize);
				PostAlignDecision(maxExtScore, reads[rid], reads[cid], b_pars, ratiophi, val->count, vss[ithread], 
					outputted, numBasesAlignedTrue, numBasesAlignedFalse, passed, matches);

				numBasesAlignedThread += endPositionV(maxExtScore.seed)-beginPositionV(maxExtScore.seed);
			}
			else // if skipAlignment == false do alignment, else save just some info on the pair to file
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
	} // all columns from start...end (omp for loop)

	double outputting = omp_get_wtime();

	int64_t* bytes = new int64_t[numThreads];
	for(int i = 0; i < numThreads; ++i)
	{
		vss[i].seekg(0, ios::end);
		bytes[i] = vss[i].tellg();
		vss[i].seekg(0, ios::beg);
	}
	int64_t bytestotal = std::accumulate(bytes, bytes+numThreads, static_cast<int64_t>(0));

	std::ofstream ofs(filename, std::ios::binary | std::ios::app);
#ifdef PRINT
	cout << "Creating output file with	" << (double)bytestotal/(double)(1024 * 1024) << " MB" << endl;
#endif
	ofs.seekp(bytestotal - 1);
	ofs.write("", 1); // this will likely create a sparse file so the actual disks won't spin yet
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
		fseek (ffinal, bytesuntil, SEEK_SET);
		std::string text = vss[ithread].str();
		fwrite(text.c_str(),1, bytes[ithread] ,ffinal);
		fflush(ffinal);
		fclose(ffinal);
	}
	delete [] bytes;
	double timeoutputt = omp_get_wtime()-outputting;

	return make_tuple(alignedpairs, alignedbases, totalreadlen, totaloutputt, totsuccbases, totfailbases, timeoutputt);
}

/**
  * Sparse multithreaded GEMM.
 **/
template <typename IT, typename NT, typename FT, typename MultiplyOperation, typename AddOperation>
void HashSpGEMM(const CSC<IT,NT>& A, const CSC<IT,NT>& B, MultiplyOperation multop, AddOperation addop, const readVector_& reads, 
	FT& getvaluetype, char* filename, const BELLApars& b_pars, const double& ratiophi)
{
	double free_memory = estimateMemory(b_pars);

#ifdef PRINT
	cout << "RAM is assumed to be:	" << free_memory / (1024 * 1024) << " MB" << endl;
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
	cout << "FLOPS is\t" << flops << endl;
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
	std::cout << "nnz (output):	"			<< nnzc					<< std::endl;
	std::cout << "Compression ratio:	"	<< compression_ratio	<< std::endl;
	std::cout << "Free memory:	"		<< free_memory		<< std::endl;
	std::cout << "Required memory:	"	<< required_memory	<< std::endl; 
	std::cout << "SpGEMM stages:	"	<< stages			<< std::endl;
//	std::cout << "max nnz per stage:	"	<< nnzcperstage		<< std::endl;
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
		double alnlenl = omp_get_wtime();
#endif
		vector<IT> * RowIdsofC = new vector<IT>[colStart[b+1]-colStart[b]];    // row ids for each column of C (bunch of cols)
		vector<FT> * ValuesofC = new vector<FT>[colStart[b+1]-colStart[b]];    // values for each column of C (bunch of cols)

		LocalSpGEMM(colStart[b], colStart[b+1], A, B, multop, addop, RowIdsofC, ValuesofC, colptrC, true);

#ifdef TIMESTEP
		double alnlen2 = omp_get_wtime();
		cout << "\nColumns [" << colStart[b] << " - " << colStart[b+1] << "] overlap time:	" << alnlen2-alnlenl << "s" << endl;
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

		tuple<size_t, size_t, size_t, size_t, size_t, size_t, double> alignstats; // (alignedpairs, alignedbases, totalreadlen, outputted, alignedtrue, alignedfalse, timeoutputt)
		alignstats = RunPairWiseAlignments(colStart[b], colStart[b+1], begnz, colptrC, rowids, values, reads, filename, b_pars, ratiophi);

#ifdef TIMESTEP
		if(!b_pars.skipAlignment)
		{
			double elapsed = omp_get_wtime()-alnlen2;
			double aligntime = elapsed-get<6>(alignstats); // substracting outputting time
			std::cout << "\nColumns ["				<< colStart[b]	<< " - "			<< colStart[b+1] << "]" << std::endl;
			std::cout << "alignmentTime:	"		<< aligntime	<< "s"				<< std::endl;
			std::cout << "alignmentRate:	"		<< (int)(static_cast<double>(get<1>(alignstats)) / aligntime) 					<< " bases/s" << std::endl;
			std::cout << "averageReadLength:	"	<< (int)(static_cast<double>(get<2>(alignstats)) / (2*get<0>(alignstats)))	<< std::endl;
			std::cout << "numPairs aligned:	"		<< get<0>(alignstats)	<< std::endl;
			cout << "averageLength of successful alignment:	"	<< (int)(static_cast<double>(get<4>(alignstats)) / get<3>(alignstats))							<< " bps" << endl;
			cout << "averageLength of failed alignment:	"		<< (int)(static_cast<double>(get<5>(alignstats)) / (get<0>(alignstats) - get<3>(alignstats)))	<< " bps" << endl;
	   }
#endif
		cout << "\nOutputted " << get<3>(alignstats) << " lines in	" << get<6>(alignstats) << "s" << endl;
		delete [] rowids;
		delete [] values;

	}//for(int b = 0; b < states; ++b)

	delete [] colptrC;
	delete [] colStart;
}

template <typename IT, typename FT>
auto RunPairWiseAlignmentsGPU(IT start, IT end, IT offset, IT * colptrC, IT * rowids, FT * values, const readVector_& reads, 
	char* filename, const BELLApars & b_pars, double ratiophi)
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

	vector<stringstream> vss(numThreads); // any chance of false sharing here? depends on how stringstream is implemented. optimize later if needed...

	vector<string> seq1s;
	vector<string> seq2s;
	vector<SeedL>  seeds;
	vector<loganResult> maxExtScoreL;

	size_t outputted = 0;

	//#pragma omp parallel for schedule(dynamic)	//	keep the order for the post evaluation code
	for(IT j = start; j < end; ++j)					//	acculate sequences for GPU batch alignment
	{
		for (IT i = colptrC[j]; i < colptrC[j+1]; ++i)
		{
			unsigned int rid = rowids[i-offset];	// row id
			unsigned int cid = j;					// column id

			const string& seq1 = reads[rid].seq;	// get reference for readibility
			const string& seq2 = reads[cid].seq;	// get reference for readibility

			unsigned short int seq1len = seq1.length();
			unsigned short int seq2len = seq2.length();

			spmatPtr_ val = values[i-offset];

			if(!b_pars.skipAlignment) // fix -z to not print 
			{
				loganResult localRes;

				//	GG: number of matching kmer into the majority voted bin
				unsigned short int matches = val->chain();
				unsigned short int minkmer = 1;

				//	GG: b_pars.minSurvivedKmers should be function of Markov chain
				if (matches < minkmer)
					continue;

				pair<int, int> kmer = val->choose();
				int i = kmer.first, j = kmer.second;

				std::string strand = "n";
				SeedL seed(i, j, i + b_pars.kmerSize, j + b_pars.kmerSize);

				std::string& seedH = seq1.substr(getBeginPositionH(seed), b_pars.kmerSize);
				std::string& seedV = seq2.substr(getBeginPositionV(seed), b_pars.kmerSize);

				std::reverse(std::begin(seedH), std::end(seedH));
				std::transform(std::begin(seedH), std::end(seedH), std::begin(seedH), complement);

				std::string cpyseq1(seq1);
				if(seedH == seedV)
				{
					strand  = "c";

					std::reverse(std::begin(cpyseq1), std::end(cpyseq1));
					std::transform(std::begin(cpyseq1), std::end(cpyseq1), std::begin(cpyseq1), complement);

					setBeginPositionH(seed, rlen - i - b_pars.kmerSize);
					setBeginPositionV(seed, j);

					setEndPositionH(seed, rlen - i);
					setEndPositionV(seed, j + b_pars.kmerSize);
				}

				loganResult localRes;
				localRes.strand = strand;

				seeds.push_back(seed);
				seq2s.push_back(seq2);
				seq1s.push_back(cpyseq1);
				maxExtScoreL.push_back(localRes);
			}
			else // if skipAlignment == false do alignment, else save just some info on the pair to file
			{
				vss[ithread] << reads[cid].nametag << '\t' << reads[rid].nametag << '\t' << val->count << '\t' << 
						seq2len << '\t' << seq1len << endl;
				++outputted;
			}
		}
	}

	if(!b_pars.skipAlignment) // fix -z to not print 
	{
		std::cout << "GPU Alignment Started" << std::endl;
		alignLogan(seq1s, seq2s, seeds, b_pars, maxExtScoreL);
		std::cout << "GPU Alignment Completed" << std::endl;

		unsigned int idx = 0;
		//	no parallelism to keep same order of pairs in alignment
		for(IT j = start; j < end; ++j) // for (end-start) columns of A^T A (one block)
		{
			size_t numAlignmentsThread   = 0;
			size_t numBasesAlignedThread = 0;
			size_t readLengthsThread     = 0;
			size_t numBasesAlignedTrue   = 0;
			size_t numBasesAlignedFalse  = 0;

			//	size_t outputted = 0;	//	moved up
			int ithread = omp_get_thread_num();

			for (IT i = colptrC[j]; i < colptrC[j+1]; ++i)	// all nonzeros in that column of A^T A
			{
				unsigned int rid = rowids[i-offset];		// row id
				unsigned int cid = j;						// column id

				const string& seq1 = reads[rid].seq;		// get reference for readibility
				const string& seq2 = reads[cid].seq;		// get reference for readibility

				unsigned short int seq1len = seq1.length();
				unsigned short int seq2len = seq2.length();

				spmatPtr_ val = values[i-offset];

#ifdef TIMESTEP
				numAlignmentsThread++;
				readLengthsThread = readLengthsThread + seq1len + seq2len;
#endif
				bool passed = false;
				loganResult maxExtScore = maxExtScoreL[idx];

				PostAlignDecisionGPU(maxExtScore, reads[rid], reads[cid], b_pars, ratiophi, val->count, 
					vss[ithread], outputted, numBasesAlignedTrue, numBasesAlignedFalse, passed);
				idx++;
#ifdef TIMESTEP
				numBasesAlignedThread += getEndPositionV(maxExtScore.seed) - getBeginPositionV(maxExtScore.seed);
#endif
			}	// all nonzeros in that column of A^T A
		}
#ifdef TIMESTEP
	//	#pragma omp critical
	//	{
			alignedpairs += numAlignmentsThread;
			alignedbases += numBasesAlignedThread;
			totalreadlen += readLengthsThread;
			totaloutputt += outputted;
			totsuccbases += numBasesAlignedTrue;
			totfailbases += numBasesAlignedFalse;
	//	}
#endif
	}	// all columns from start...end (omp for loop)

	double outputting = omp_get_wtime();

	int64_t * bytes = new int64_t[numThreads];
	for(int i=0; i< numThreads; ++i)
	{
		vss[i].seekg(0, ios::end);
		bytes[i] = vss[i].tellg();
		vss[i].seekg(0, ios::beg);
	}
	int64_t bytestotal = std::accumulate(bytes, bytes+numThreads, static_cast<int64_t>(0));

	std::ofstream ofs(filename, std::ios::binary | std::ios::app);
#ifdef PRINT
	cout << "Creating or appending to output file with " << (double)bytestotal/(double)(1024 * 1024) << " MB" << endl;
#endif
	ofs.seekp(bytestotal - 1);
	ofs.write("", 1); // this will likely create a sparse file so the actual disks won't spin yet
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
	double timeoutputt = omp_get_wtime()-outputting;
	return make_tuple(alignedpairs, alignedbases, totalreadlen, totaloutputt, totsuccbases, totfailbases, timeoutputt);
}

/**
  * Sparse multithreaded GEMM.
 **/
template <typename IT, typename NT, typename FT, typename MultiplyOperation, typename AddOperation>
void HashSpGEMMGPU(const CSC<IT,NT> & A, const CSC<IT,NT> & B, MultiplyOperation multop, AddOperation addop, const readVector_ & reads, 
	FT & getvaluetype, char* filename, const BELLApars & b_pars, double ratiophi)
{
	double free_memory = estimateMemory(b_pars);

#ifdef PRINT
	cout << "RAM is assumed to be:	" << free_memory / (1024 * 1024) << " MB" << endl;
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
	cout << "FLOPS is\t" << flops << endl;
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
	std::cout << "nnz (output):	"			<< nnzc					<< std::endl;
	std::cout << "Compression ratio:	"	<< compression_ratio	<< std::endl;
	std::cout << "Free memory:	"		<< free_memory		<< std::endl;
	std::cout << "Required memory:	"	<< required_memory	<< std::endl; 
	std::cout << "SpGEMM stages:	"	<< stages			<< std::endl;
	//	std::cout << "max nnz per stage:	"	<< nnzcperstage		<< std::endl;
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
		double alnlen2 = omp_get_wtime();
		cout << "\nColumns [" << colStart[b] << " - " << colStart[b+1] << "] overlap time:	" << alnlen2-alnlenl << "s" << endl;
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

		tuple<size_t, size_t, size_t, size_t, size_t, size_t, double> alignstats; // (alignedpairs, alignedbases, totalreadlen, outputted, alignedtrue, alignedfalse, timeoutputt)
		alignstats = RunPairWiseAlignmentsGPU(colStart[b], colStart[b+1], begnz, colptrC, rowids, values, reads, b_pars, filename, ratiophi);

#ifdef TIMESTEP
		if(!b_pars.skipAlignment)
		{
			double elapsed = omp_get_wtime()-alnlen2;
			double aligntime = elapsed-get<6>(alignstats); // substracting outputting time
			std::cout << "\nColumns ["				<< colStart[b]	<< " - "			<< colStart[b+1] << "]" << std::endl;
			std::cout << "alignmentTime:	"		<< aligntime	<< "s"				<< std::endl;
			std::cout << "alignmentRate:	"		<< (int)(static_cast<double>(get<1>(alignstats)) / aligntime) 					<< " bases/s" << std::endl;
			std::cout << "averageReadLength:	"	<< (int)(static_cast<double>(get<2>(alignstats)) / (2*get<0>(alignstats)))	<< std::endl;
			std::cout << "numPairs aligned:	"		<< get<0>(alignstats)	<< std::endl;
			cout << "averageLength of successful alignment:	"	<< (int)(static_cast<double>(get<4>(alignstats)) / get<3>(alignstats))							<< " bps" << endl;
			cout << "averageLength of failed alignment:	"		<< (int)(static_cast<double>(get<5>(alignstats)) / (get<0>(alignstats) - get<3>(alignstats)))	<< " bps" << endl;
		}
#endif
		cout << "\nOutputted " << get<3>(alignstats) << " lines in " << get<6>(alignstats) << "s" << endl;
		delete [] rowids;
		delete [] values;

	}//for(int b = 0; b < states; ++b)

	delete [] colptrC;
	delete [] colStart;
}
