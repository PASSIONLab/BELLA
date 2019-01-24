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
	bool adapThr; 			// Apply adaptive alignment threshold (v)
	int defaultThr;   		// default alignment score threshold (a), only matters when adapThr=false, to be deprecated	
	bool alignEnd;			// Filter out alignments not achieving end of the read "relaxed" (x)
	int relaxMargin;		// epsilon parameter for alignment on edges (w)
	double deltaChernoff;	// delta computed via Chernoff bound (c)

	BELLApars():totalMemory(8000.0), userDefMem(false), kmerRift(1000), skipEstimate(false), skipAlignment(false), adapThr(true), defaultThr(50),
			alignEnd(false), relaxMargin(300), deltaChernoff(0.1) {};
};

template <typename T>
    bool isinrift(const T& value, const T& left, const T& right) {
    return (value > left) && (value < right);
}

typedef seqan::Seed<seqan::Simple> TSeed;
struct seqAnResult {
    int score;
    string strand;
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

    int count = 0;   // number of shared k-mers
    int pos[4] = {0};  // pos1i, pos1j, pos2i, pos2j 
};

typedef shared_ptr<spmatType_> spmatPtr_; // pointer to spmatType_ datastruct
typedef std::vector<Kmer> Kmers;

struct alignmentInfo {

    int64_t score;              // score 
    uint32_t apos, bpos;        // (8) pos in the sections 
    uint32_t alen, blen;        // (8) lengths of the segments 
};

#endif