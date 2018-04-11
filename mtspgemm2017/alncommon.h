#ifndef ALN_COMMON_H_
#define ALN_COMMON_H_

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/seeds.h>
using namespace seqan;

typedef Seed<Simple> TSeed;
struct seqAnResult {
    int score;
    string strand;
    TSeed seed;
};

#endif