#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/seeds.h>
#include "common.h"
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
#include "../kmercode/Kmer.hpp"

using namespace seqan;
using namespace std;

double adaptiveSlope(double error)
{
    double p_mat = pow(1-error,2);  // match
    double p_mis = 1-p_mat;         // mismatch/gap
    double alpha = 1;               // match penalty
    double beta = 1;                // mismatch/gap penalty

    return alpha*p_mat - beta*p_mis;
}

bool toEnd(int colStart, int colEnd, int colLen, int rowStart, int rowEnd, int rowLen, int relaxMargin)
{
    int minLeft = min(colStart, rowStart);
    int minRight = min(colLen-colEnd, rowLen-rowEnd);

     if(minLeft-relaxMargin <= 0)
        minLeft = 0;
     if(minRight-relaxMargin <= 0)
        minRight = 0;

     if((minLeft == 0 || minRight == 0))
        return true;
    else
        return false;
}

/**
 * @brief alignSeqAn does the seed-and-extend alignment
 * @param row
 * @param col
 * @param rlen is the length of the row sequence
 * @param i is the starting position of the k-mer on the first read
 * @param j is the starting position of the k-mer on the second read
 * @param xdrop
 * @param useHOPC specifies whether HOPC representations of kmers are used
 * @return alignment score and extended seed
 */
seqAnResult alignSeqAn(const std::string & row, const std::string & col, int rlen, int i, int j, int xdrop, int kmer_len, bool useHOPC, bool iRev, bool jRev) {

    Score<int, Simple> scoringScheme(1,-1,-1);

    Dna5String seqH(row);
    Dna5String seqV(col);
    Dna5String seedH;
    Dna5String seedV;
    string strand;
    int longestExtensionTemp;
    seqAnResult longestExtensionScore;


    TSeed seed(i, j, i+kmer_len, j+kmer_len);
    seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
    seedV = infix(seqV, beginPositionV(seed), endPositionV(seed));

    /* we are reversing the "row", "col" is always on the forward strand */
    Dna5StringReverseComplement twin(seedH);

    bool hopcReverse = useHOPC && ( iRev != jRev );

    if ( (twin == seedV) || hopcReverse )
    {
        strand = 'c';
        Dna5StringReverseComplement twinRead(seqH);
        i = rlen-i-kmer_len;

        setBeginPositionH(seed, i);
        setBeginPositionV(seed, j);
        setEndPositionH(seed, i+kmer_len);
        setEndPositionV(seed, j+kmer_len);

        /* Perform match extension */
        longestExtensionTemp = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, xdrop, kmer_len, GappedXDrop());

    } else
    {
        strand = 'n';
        longestExtensionTemp = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, xdrop, kmer_len, GappedXDrop());
    }

    longestExtensionScore.score = longestExtensionTemp;
    longestExtensionScore.seed = seed;
    longestExtensionScore.strand = strand;
    return longestExtensionScore;
}

#endif
