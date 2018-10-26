#ifndef _GLOBAL_H_
#define _GLOBAL_H_

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

    int count = 0;   /* number of shared k-mers */
    int pos[4] = {0};  /* pos1i, pos1j, pos2i, pos2j */
};

typedef shared_ptr<spmatType_> spmatPtr_; // pointer to spmatType_ datastruct
typedef std::vector<Kmer> Kmers;

struct alignmentInfo {

    int64_t score;              /** score */
    uint32_t apos, bpos;        /** (8) pos in the sections */
    uint32_t alen, blen;        /** (8) lengths of the segments */
};

#ifdef __cplusplus
extern "C" {
#endif

#include "../optlist/optlist.h" // command line parser
#include "gaba.h" // sequence alignment lib

#ifdef __cplusplus
}
#endif

#include "../libcuckoo/cuckoohash_map.hh"

#endif