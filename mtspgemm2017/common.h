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

    int kmerRift;
	bool skipEstimate;  	// Do not estimate error but use user-defined error (e)
	bool skipAlignment;  	// Do not align (z)
    bool allKmer;           // Use all possible kmers (non-overlapping and separated by <kmerRift> bases) as alignment seeds (K)
	bool adapThr; 			// Apply adaptive alignment threshold (v)
	int defaultThr;   		// default alignment score threshold (a), only matters when adapThr=false, to be deprecated	
	bool alignEnd;			// Filter out alignments not achieving end of the read "relaxed" (x)
	int relaxMargin;		// epsilon parameter for alignment on edges (w)
	double deltaChernoff;	// delta computed via Chernoff bound (c)
    bool outputPaf;         // output in paf format (p)

	BELLApars():totalMemory(8000.0), userDefMem(false), kmerRift(1000), skipEstimate(false), skipAlignment(false), allKmer(false), adapThr(true), defaultThr(50),
			alignEnd(false), relaxMargin(300), deltaChernoff(0.2), outputPaf(false) {};
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

struct spmatType_ {

    int count = 0;              // number of shared k-mers
    vector<pair<int,int>> pos;  // vector of k-mer positions <read-i, read-j> (if !K, use at most 2 kmers, otherwise all)
};

typedef shared_ptr<spmatType_> spmatPtr_; // pointer to spmatType_ datastruct

typedef std::vector<Kmer> Kmers;

struct alignmentInfo {

    int64_t score;              // score 
    uint32_t apos, bpos;        // (8) pos in the sections 
    uint32_t alen, blen;        // (8) lengths of the segments 
};

#endif