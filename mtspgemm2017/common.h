#ifndef _COMMON_H_
#define _COMMON_H_

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/seeds.h>
#ifndef __NVCC__
	#include "../xavier/xavier.h"
#endif

#ifndef PRINT
#define PRINT
#endif

#include "../kmercode/hash_funcs.h"
#include "../kmercode/Kmer.hpp"
#include "../kmercode/Buffer.h"
#include "../kmercode/common.h"
#include "../kmercode/fq_reader.h"
#include "../kmercode/ParallelFASTQ.h"
#include "../kmercode/bound.hpp"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

#include "../optlist/optlist.h" // command line parser

#ifdef __cplusplus
}
#endif

#include "../libcuckoo/cuckoohash_map.hh"

#ifdef PRINT
	#define printLog(var) do { std::cerr << "INFO:	" << __FILE__ << "(" << __LINE__ << ")	" << #var << " = " << (var) << std::endl; } while(0)
#else
	#define printLog(var)
#endif

struct BELLApars
{
	unsigned short int		kmerSize;			// KmerSize
	unsigned short int		binSize;			// Bin size chaining algorithm 			(w)
	unsigned short int		minOverlap;			// minimum overlap length	 			(l)
	short int		        fixedThreshold;		// Default alignment score threshold 	(a)
	unsigned short int		xDrop;				// SeqAn xDrop value 					(x)
	unsigned short int		numGPU;				// Number GPUs available/to be used  	(g)
	short int		        myMarkovOverlap;	// Overlap length to see 1 shared k-mer (b)
	bool	skipEstimate;		// Do not estimate error but use user-defined error 	(e)
	bool	skipAlignment;		// Do not align 										(z)
	bool	outputPaf;			// Output in paf format 								(p)
	bool	userDefMem;			// RAM available 										(m)
	double	deltaChernoff;		// delta computed via Chernoff bound 					(d)
	double	totalMemory;		// In MB, default is ~ 8GB
	double	errorRate;			// default error rate if estimation is disable 			(e)
	double	minProbability;		// reliable range probability threshold 				(r)
	double	minpNMC;			// nested markov chain probability threshold 		    (n)

	BELLApars(): kmerSize(17), binSize(500), minOverlap(2000), fixedThreshold(-1), xDrop(7), numGPU(1), 
					myMarkovOverlap(-1), skipEstimate(true), skipAlignment(false), outputPaf(false), userDefMem(false),
						deltaChernoff(0.10), totalMemory(8000.0), errorRate(0.00), minProbability(0.10), minpNMC(0.90) {};
};

template <typename T>
	bool isinrift(const T& value, const T& left, const T& right) {
	return (value > left) && (value < right);
}

#ifndef __NVCC__

struct xavierResult {
    int score;
    std::string strand;
    SeedX seed;
};

#endif

typedef seqan::Seed<seqan::Simple> TSeed;
struct seqAnResult {
	int score;
	std::string strand;
	TSeed seed;
};

struct readType_ {
	std::string nametag;
	std::string seq; 
	int readid;

	bool operator < (readType_ & str)
	{
		return (readid < str.readid);
	}
};

typedef std::vector<readType_> readVector_;

// EK: sort function for sorting a std::vector of indices by the values in a std::vector of int
struct SortBy:std::binary_function<unsigned short int, unsigned short int, bool>
{
	SortBy(const std::vector<unsigned short int>& par) : vec(par) {}
	bool operator()(unsigned short int idx1, unsigned short int idx2) const { return vec[idx1] > vec[idx2]; }
	const std::vector<unsigned short int>& vec;
};

struct spmatType_ {

	unsigned short int count = 0;		// number of shared k-mers
	std::vector<std::vector<pair<unsigned short int, unsigned short int>>> pos;	// std::vector of k-mer positions <read-i, read-j> (if !K, use at most 2 kmers, otherwise all) per bin
	std::vector<unsigned short int> support;	// number of k-mers supporting a given overlap
	std::vector<unsigned short int> overlap;	// overlap values
	std::vector<unsigned short int> ids;		// indices corresponded to sorting of support

	//	GG: sort according to support number
	void sort() {
		ids = std::vector<unsigned short int>(support.size());					// number of support
		std::iota(ids.begin(), ids.end(), 0);				// assign an id
		std::sort(ids.begin(), ids.end(), SortBy(support));	// sort support by supporting k-mers
	}

	//	GG: print overlap estimate and support number
	void print() {
		std::copy(overlap.begin(), overlap.end(), std::ostream_iterator<unsigned short int>(std::cout, "\t")); std::cout << std::endl;
		std::copy(support.begin(), support.end(), std::ostream_iterator<unsigned short int>(std::cout, "\t")); std::cout << std::endl;
	}

	//	GG: number of kmer supporting the most voted bin
	int chain() {
		ids = std::vector<unsigned short int>(support.size());	// number of support
		std::iota(ids.begin(), ids.end(), 0);				// assign an id
		std::sort(ids.begin(), ids.end(), SortBy(support));	// sort support by supporting k-mers

		ids.resize(1);				// GG: we don't care about other support, we want only the majority voted one
		return support[ids[0]];		// number of kmer in the most voted bin
	}

	//	GG: overlap len in the most voted bin
	int overlaplength() {
		ids = std::vector<unsigned short int>(support.size());	// number of support
		std::iota(ids.begin(), ids.end(), 0);				// assign an id
		std::sort(ids.begin(), ids.end(), SortBy(support));	// sort support by supporting k-mers

		ids.resize(1);				// GG: we don't care about other support, we want only the majority voted one
		return overlap[ids[0]];		// number of kmer in the most voted bin
	}

	//	GG: choose does also sorting and return the position of the first k-mer in each bin
	pair<unsigned short int, unsigned short int> choose() {
		ids = std::vector<unsigned short int>(support.size());	// number of support
		std::iota(ids.begin(), ids.end(), 0);				// assign an id
		std::sort(ids.begin(), ids.end(), SortBy(support));	// sort support by supporting k-mers

		ids.resize(1);			// GG: we don't care about other support, we want only the majority voted one
		pos[ids[0]].resize(1);	// GG: same for the number of kmers in the choosen bin, we need only one

		return pos[ids[0]][0];	// GG: returning choosen seed // It might be better choose the kmer randomly and find an x-drop/binsize ratio to justify
	}

};

typedef shared_ptr<spmatType_> spmatPtr_; // pointer to spmatType_ datastruct
typedef std::vector<Kmer> Kmers;

struct alignmentInfo {

	int score;	//	score
	unsigned short int apos, bpos;	// (8)	pos in the sections
	unsigned short int alen, blen;	// (8)	lengths of the segments
};

// modified from https://www.geeksforgeeks.org/interval-tree/
struct Interval 
{ 
    int lower, upper, num;
}; 
  
// node struct in interval search tree 
struct ITNode 
{
    Interval *i;
    ITNode *le, *ri; 
}; 
  
// function to create a new interval search tree node 
ITNode * create(Interval i) 
{ 
	ITNode *temp = new ITNode; 
	temp->i = new Interval(i);
	temp->le = temp->ri = NULL; 
}; 

// function to insert a new interval search tree node 
ITNode *insert(ITNode *root, Interval i) 
{ 
    // tree is empty, new node becomes root 
    if (root == NULL) 
        return create(i); 
  
    // low value of interval at root 
    int l = root->i->lower; 
  
    if (i.lower < l) root->le = insert(root->le, i);
    else root->ri = insert(root->ri, i);
  
    return root; 
}

// function to check if L is in overlap
bool within(Interval i, int L) 
{ 
    if (i.lower <= L && L < i.upper) return true; 
    return false; 
}

// main function that searches a given interval i in a given interval tree 
Interval *search(ITNode *root, int L) 
{ 
    // Base case the tree is empty
    if (root == NULL) return NULL; 
  
	// L is in this interval
    if (within(*(root->i), L)) 
		return root->i; 

	// L is greater and we are in the ri-most node
	if (root->ri == NULL && L >= root->i.upper) 
	{
		root->i.num++;
		return root->i;	
	}

	// L is smaller and we are can go left
    if (root->le != NULL && L < root->i.lower) 
		return search(root->le, L);
	
	// L is greater and we are can go left
	if (root->ri != NULL && L >= root->i.upper) 
		return search(root->ri, L); 
}

#endif
