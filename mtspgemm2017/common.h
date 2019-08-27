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
	unsigned short int		kmerSize;			// KmerSize
	unsigned short int		kmerRift;			// minDistance between Kmer
	unsigned short int		minOverlap;			// minOverlap length to detect (used to select the number number of shared kmer)
	unsigned short int		minSurvivedKmers;	// GG: to be mathematically determine via Markov chain with minOverlap and error rate
	unsigned short int		maxOverhang;		// maxOverhang
	unsigned short int		maxJump;			// maxJump to detect chimeric sequences
	unsigned short int		binSize;			// bin size chaining algorithm (b)
	unsigned short int		defaultThr;			// default alignment score threshold (a), only matters when adapThr=false, GG: to be deprecated
	unsigned short int		xDrop;				// seqAn xDrop value (7)
	bool	skipEstimate;		// Do not estimate error but use user-defined error (e)
	bool	skipAlignment;		// Do not align (z)
	bool	adapThr;			// Apply adaptive alignment threshold (v)
	bool	outputPaf;			// output in paf format (p)
	bool	userDefMem;
	bool	useGerbil;
	bool	enableGPU;
	float	maxDivergence;		// maxDivergence to output a pair
	double	deltaChernoff;		// delta computed via Chernoff bound (c)
	double	totalMemory;		// in MB, default is ~ 8GB
	double	errorRate;			// default error rate if estimation is disable (e)
	double	minProbability;		// reliable range probability threshold (y)

	BELLApars(): kmerSize(17), kmerRift(kmerSize), minOverlap(1000), minSurvivedKmers(-1), maxOverhang(1500), maxJump(1500), binSize(500), defaultThr(0),
					xDrop(7), skipEstimate(false), skipAlignment(false), adapThr(true), outputPaf(false), userDefMem(false), useGerbil(false), enableGPU(false), maxDivergence(0.25),
						deltaChernoff(0.10), totalMemory(8000.0), errorRate(0.15) {};
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
struct SortBy:std::binary_function<unsigned short int, unsigned short int, bool>
{
	SortBy(const std::vector<unsigned short int>& par) : vec(par) {}
	bool operator()(unsigned short int idx1, unsigned short int idx2) const { return vec[idx1] > vec[idx2]; }
	const std::vector<unsigned short int>& vec;
};

struct spmatType_ {

	unsigned short int count = 0;		// number of shared k-mers
	vector<vector<pair<unsigned short int, unsigned short int>>> pos;	// vector of k-mer positions <read-i, read-j> (if !K, use at most 2 kmers, otherwise all) per bin
	vector<unsigned short int> support;	// number of k-mers supporting a given overlap
	vector<unsigned short int> overlap;	// overlap values
	vector<unsigned short int> ids;		// indices corresponded to sorting of support

	//	GG: sort according to support number
	void sort() {
		ids = vector<unsigned short int>(support.size());					// number of support
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
		ids = vector<unsigned short int>(support.size());	// number of support
		std::iota(ids.begin(), ids.end(), 0);				// assign an id
		std::sort(ids.begin(), ids.end(), SortBy(support));	// sort support by supporting k-mers

		ids.resize(1);				// GG: we don't care about other support, we want only the majority voted one
		return support[ids[0]];		// number of kmer in the most voted bin
	}

	//	GG: overlap len in the most voted bin
	int overlaplength() {
		ids = vector<unsigned short int>(support.size());	// number of support
		std::iota(ids.begin(), ids.end(), 0);				// assign an id
		std::sort(ids.begin(), ids.end(), SortBy(support));	// sort support by supporting k-mers

		ids.resize(1);				// GG: we don't care about other support, we want only the majority voted one
		return overlap[ids[0]];		// number of kmer in the most voted bin
	}

	//	GG: choose does also sorting and return the position of the first k-mer in each bin
	pair<unsigned short int, unsigned short int> choose() {
		ids = vector<unsigned short int>(support.size());	// number of support
		std::iota(ids.begin(), ids.end(), 0);				// assign an id
		std::sort(ids.begin(), ids.end(), SortBy(support));	// sort support by supporting k-mers

		ids.resize(1);			// GG: we don't care about other support, we want only the majority voted one
		pos[ids[0]].resize(1);	// GG: same for the number of kmers in the choosen bin, we need only one

		return pos[ids[0]][0];	// GG: returning choosen seed
	}

};

typedef shared_ptr<spmatType_> spmatPtr_; // pointer to spmatType_ datastruct
typedef std::vector<Kmer> Kmers;

struct alignmentInfo {

	int score;	//	score
	unsigned short int apos, bpos;	// (8)	pos in the sections
	unsigned short int alen, blen;	// (8)	lengths of the segments
};

#endif
