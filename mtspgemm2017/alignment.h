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
// #define _UNGAPPED

typedef Seed<Simple>  TSeed;
typedef SeedSet<TSeed> TSeedSet;

pair<int64_t,Seed<Simple>> seqanAlOne(std::string & row, std::string & col, int rlen, int i, int j, int dropFactor) {

	Score<int, Simple> scoringScheme(1, -1, -1);

	Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    int64_t longestExtensionTemp;
    std::pair<int64_t,Seed<Simple>> longestExtensionScore;

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
        longestExtensionTemp = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        #endif

        #ifdef _UNGAPPED
        longestExtensionTemp = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        #endif
    } else
    {
        /* Perform match extension */
        #ifdef _GAPPED
        longestExtensionTemp = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        #endif

        #ifdef _UNGAPPED
        longestExtensionTemp = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        #endif
    }

    longestExtensionScore = make_pair(longestExtensionTemp, seed1);
    return longestExtensionScore;
}

pair<int64_t,Seed<Simple>> seqanAlGen(std::string & row, std::string & col, int rlen, int i, int j, int l, int m, int dropFactor) {

	Score<int, Simple> scoringScheme(1, -1, -1);

	Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    std::pair<int64_t, int64_t> longestExtensionTemp;
    std::pair<int64_t,Seed<Simple>> longestExtensionScore;

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
        longestExtensionTemp.first = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        longestExtensionTemp.second = extendSeed(seed2, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        #endif

        #ifdef _UNGAPPED
        longestExtensionTemp.first = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        longestExtensionTemp.second = extendSeed(seed2, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        #endif
    } else
    {
        #ifdef _GAPPED
        longestExtensionTemp.first = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        longestExtensionTemp.second = extendSeed(seed2, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        #endif

        #ifdef _UNGAPPED
        longestExtensionTemp.first = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        longestExtensionTemp.second = extendSeed(seed2, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        #endif
    }

    if(longestExtensionTemp.first > longestExtensionTemp.second)
    {
        longestExtensionScore = make_pair(longestExtensionTemp.first, seed1);
    } else
    {
        longestExtensionScore = make_pair(longestExtensionTemp.second, seed2);
    }

    return longestExtensionScore;
}

pair<int64_t,Seed<Simple>> seqanAlOneAllKmer(std::string & row, std::string & col, int rlen, std::vector<std::pair<int,int>> vpos, int dropFactor) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    int64_t longestExtensionTemp;
    std::pair<int64_t,Seed<Simple>> longestExtensionScore;

    seqH = row;
    seqV = col;

    std::vector<std::pair<int,int>>::iterator it;
    it = vpos.begin();

    Seed<Simple> seed(it->first, it->second, it->first+KMER_LENGTH, it->second+KMER_LENGTH);
    seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
    seedV = infix(seqV, beginPositionV(seed), endPositionV(seed));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        Dna5StringReverseComplement twinRead(seqH);
        int i = it->first;
        i = rlen-i-KMER_LENGTH;
        Seed<Simple> seed(i, it->second, i+KMER_LENGTH, it->second+KMER_LENGTH);

        /* Perform match extension */
        #ifdef _GAPPED
        longestExtensionTemp = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        #endif

        #ifdef _UNGAPPED
        longestExtensionTemp = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        #endif
    } else
    {
        /* Perform match extension */
        #ifdef _GAPPED
        longestExtensionTemp = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        #endif

        #ifdef _UNGAPPED
        longestExtensionTemp = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
        #endif
    }

    longestExtensionScore = make_pair(longestExtensionTemp, seed);
    return longestExtensionScore;
}

pair<int64_t,Seed<Simple>> seqanAlGenAllKmer(std::string & row, std::string & col, int rlen, std::vector<std::pair<int,int>> vpos, int dropFactor) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    std::pair<int64_t,Seed<Simple>> tempScore;
    std::pair<int64_t,Seed<Simple>> longestExtensionScore;

    seqH = row;
    seqV = col;
    longestExtensionScore.first = 0;

    for(std::vector<std::pair<int,int>>::iterator it=vpos.begin(); it!=vpos.end(); ++it)
    {
        Seed<Simple> seed(it->first, it->second, it->first+KMER_LENGTH, it->second+KMER_LENGTH);
        seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
    
        Dna5StringReverseComplement twin(seedH);
    
        if(twin == seedV)
        {
            Dna5StringReverseComplement twinRead(seqH);

            int i = it->first;
            i = rlen-i-KMER_LENGTH;
            Seed<Simple> seed(i, it->second, i+KMER_LENGTH, it->second+KMER_LENGTH);
        
            /* Perform match extension */
            #ifdef _GAPPED
            tempScore.first = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
            #endif
    
            #ifdef _UNGAPPED
            tempScore.first = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
            #endif
        } else
        {
            #ifdef _GAPPED
            tempScore.first = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
            #endif
    
            #ifdef _UNGAPPED
            tempScore.first = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, UnGappedXDrop());
            #endif
        }
        
        if(tempScore.first > longestExtensionScore.first)
            longestExtensionScore = make_pair(tempScore.first, seed);
    }
    return longestExtensionScore;
}