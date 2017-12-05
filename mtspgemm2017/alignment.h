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
using namespace std;

#define KMER_LENGTH 17
#define _GAPPED
//#define _UNGAPPED

typedef Seed<Simple>  TSeed;
typedef SeedSet<TSeed> TSeedSet;

int64_t seqanAlOne(std::string & row, std::string & col, int rlen, int i, int j, int dropFactor) {

	Score<int, Simple> scoringScheme(1, -1, -1);

	Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    int64_t longestExtensionScoreOne;

	seqH = row;
    seqV = col;

    Seed<Simple> seed1(i, j, i+KMER_LENGTH, j+KMER_LENGTH);
    seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
    seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        Dna5StringReverseComplement twinRead(seqH);
        i = rlen-i-KMER_LENGTH;
        Seed<Simple> seed1(i, j, i+KMER_LENGTH, j+KMER_LENGTH);

        /* Perform match extension */
        #ifdef _GAPPED
        longestExtensionScoreOne = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        #endif

        #ifdef _UNGAPPED
        longestExtensionScoreOne = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        #endif
    } else
    {
        /* Perform match extension */
        #ifdef _GAPPED
        longestExtensionScoreOne = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        #endif

        #ifdef _UNGAPPED
        longestExtensionScoreOne = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        #endif
    }

    return longestExtensionScoreOne;
}

pair<int64_t, int64_t> seqanAlGen(std::string & row, std::string & col, int rlen, int i, int j, int l, int m, int dropFactor) {

	Score<int, Simple> scoringScheme(1, -1, -1);

	Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    pair<int64_t, int64_t> longestExtensionScore;

	seqH = row;
    seqV = col;

    Seed<Simple> seed1(i, j, i+KMER_LENGTH, j+KMER_LENGTH);
    Seed<Simple> seed2(l, m, l+KMER_LENGTH, m+KMER_LENGTH);
    seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
    seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        Dna5StringReverseComplement twinRead(seqH);

        i = rlen-i-KMER_LENGTH;
        l = rlen-l-KMER_LENGTH;

        Seed<Simple> seed1(i, j, i+KMER_LENGTH, j+KMER_LENGTH);
        Seed<Simple> seed2(l, m, l+KMER_LENGTH, m+KMER_LENGTH);

        /* Perform match extension */
        #ifdef _GAPPED
        longestExtensionScore.first = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        longestExtensionScore.second = extendSeed(seed2, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        #endif

        #ifdef _UNGAPPED
        longestExtensionScore.first = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        longestExtensionScore.second = extendSeed(seed2, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        #endif
    } else
    {
        #ifdef _GAPPED
        longestExtensionScore.first = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        longestExtensionScore.second = extendSeed(seed2, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        #endif

        #ifdef _UNGAPPED
        longestExtensionScore.first = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        longestExtensionScore.second = extendSeed(seed2, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        #endif
    }
    return longestExtensionScore;
}