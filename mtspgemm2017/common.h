#ifndef _COMMON_H_
#define _COMMON_H_

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/seeds.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "../optlist/optlist.h" // command line parser
//#include "gaba.h" // sequence alignment libgaba

#ifdef __cplusplus
}
#endif

#include "../libcuckoo/cuckoohash_map.hh"

struct BELLApars
{
	double totalMemory;	// in MB, default is ~ 8GB
	bool userDefMem;

	int  kmerRift;
	bool skipEstimate;  	// Do not estimate error but use user-defined error (e)
	bool skipAlignment;  	// Do not align (z)
	bool allKmer;           // Use all possible kmers (non-overlapping and separated by <kmerRift> bases) as alignment seeds (K)
	bool adapThr; 			// Apply adaptive alignment threshold (v)
	int  defaultThr;   		// default alignment score threshold (a), only matters when adapThr=false, to be deprecated	
	bool alignEnd;			// Filter out alignments not achieving end of the read "relaxed" (x)
	int  relaxMargin;		// epsilon parameter for alignment on edges (w)
	bool outputPaf;         // output in paf format (p)
	int  bin;				// bin size chaining algorithm (b)
	double deltaChernoff;	// delta computed via Chernoff bound (c)

	BELLApars():totalMemory(8000.0), userDefMem(false), kmerRift(1000), skipEstimate(false), skipAlignment(false), allKmer(false), adapThr(true), defaultThr(50),
			alignEnd(false), relaxMargin(300), outputPaf(false), bin(500), deltaChernoff(0.2) {};
};

template <typename T>
	bool isinrift(const T& value, const T& left, const T& right) {
	return (value > left) && (value < right);
}

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

typedef vector<readType_> readVector_;

// EK: sort function for sorting a vector of indices by the values in a vector of int
struct SortBy:std::binary_function<int, int, bool>
{
	SortBy(const std::vector<int>& par) : vec(par) {}
	bool operator()(int idx1, int idx2) const { return vec[idx1] > vec[idx2]; }
	const std::vector<int>& vec;
};

struct spmatType_ {

	int count = 0;						// number of shared k-mers
	vector<vector<pair<int, int>>> pos;	// vector of k-mer positions <read-i, read-j> (if !K, use at most 2 kmers, otherwise all) per bin
	vector<int> support;				// number of k-mers supporting a given overlap
	vector<int> overlap; 				// overlap values
	vector<int> ids;					// indices corresponded to sorting of support (GG:?)

	//	GG: debug
	void sort() {
		ids = vector<int>(support.size());					// number of support
		std::iota(ids.begin(), ids.end(), 0);				// assign an id
		std::sort(ids.begin(), ids.end(), SortBy(support));	// sort support by supporting k-mers
	}

	//	GG: debug
	void print() {
		std::cout << "OVERLAPS\t";
		std::copy(overlap.begin(), overlap.end(), std::ostream_iterator<int>(std::cout, "\t")); std::cout << std::endl;
		std::cout << "SUPPORTS\t";
		std::copy(support.begin(), support.end(), std::ostream_iterator<int>(std::cout, "\t")); std::cout << std::endl;
	}

	//	GG: choose does also sorting and return the position of the first k-mer in each bin
	//	vector<pair<int, int>> choose() {
	//	GG: choose does also sorting and return the positions all of all the k-mers
	vector<vector<pair<int, int>>> choose() {

		ids = vector<int>(support.size());					// number of support
		std::iota(ids.begin(), ids.end(), 0);				// assign an id
		std::sort(ids.begin(), ids.end(), SortBy(support));	// sort support by supporting k-mers

		//	ids.resize(1);	// GG: we don't care about other support, we want only the majority voted one
		//	vector<pair<int, int>> kmervect;
		vector<vector<pair<int, int>>> kmervect;
		for(auto i = ids.begin(); i != ids.end(); i++)	// GG: keep saved one k-mer per bin in case of bins with ties or similar number of support
			kmervect.push_back(pos[*i]);				// GG: keep saving all k-mers positions for distirbutions purpose
			//	kmervect.push_back(pos[*i][0]);			// GG: same for the number of kmers in the choosen bin, we need only one

		//	return pos[ids[0]][0];	// GG: returning choosen seed
		return kmervect;			// GG: return vector of choosen seeds of size nbins to manage ties
	}
};

typedef shared_ptr<spmatType_> spmatPtr_; // pointer to spmatType_ datastruct

typedef std::vector<Kmer> Kmers;

struct alignmentInfo {

	int64_t score;              // score 
	uint32_t apos, bpos;        // (8) pos in the sections 
	uint32_t alen, blen;        // (8) lengths of the segments 
};

#endif
