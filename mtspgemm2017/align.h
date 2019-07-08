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

bool matchCaps(int cap, bool reverse, const std::string seq1, const std::string seq2, int start1, int start2, bool useHOPC, int kmer_len) {
    // First version of cap implementation (item 24 in doc)
    // Assess whether to align the pair based on the caps on the end of the kmer
    // std::string beg1, beg2, end1, end2;
    //
    // if ( start1 < cap || start2 < cap ) {
    //     beg1 = "";
    //     beg2 = "";
    // } else {
    //     beg1 = seq1.substr(start1-cap, cap);
    //     beg2 = seq2.substr(start2-cap, cap);
    // }
    //
    // if ( ( seq1.length() - start1 - kmer_len < cap ) || ( seq2.length() - start2 - kmer_len < cap ) ) {
    //     end1 = "";
    //     end2 = "";
    // } else {
    //     end1 = seq1.substr(start1+kmer_len, cap);
    //     end2 = seq2.substr(start2+kmer_len, cap);
    // }
    //
    // if ( reverse ) {
    //     beg2.swap(end2);
    //     std::reverse(beg2.begin(), beg2.end());
    //     std::reverse(end2.begin(), end2.end());
    // }
    //
    // return (beg1 == beg2) && (end1 == end2);

    // Second version of cap implementation (item 25 in doc)
    int difs = 0;
    if ( start1 > cap && start2 > cap && ( seq1.length() - start1 - kmer_len > cap ) && ( seq2.length() - start2 - kmer_len > cap ) ) {
      std::string cap1 = seq1.substr(start1-cap, 2*cap+kmer_len);
      std::string cap2 = seq2.substr(start2-cap, 2*cap+kmer_len);
      if ( reverse ) std::reverse(cap2.begin(), cap2.end());
      for(int idx = 0; idx < 2*cap+kmer_len; idx++)
        if (cap1[idx] != cap2[idx]) difs++;
      // if ( difs > ( 2*cap ) || (!useHOPC && difs > cap ) ) {
      //   longestExtensionScore.score = 0;
      //   longestExtensionScore.seed = seed;
      //   longestExtensionScore.strand = '0';
      //   return longestExtensionScore;
      // }
      return difs < (2*cap) || (!useHOPC && difs < cap);
    }
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
seqAnResult alignSeqAn(const std::string & row, const std::string & col, int rlen, int i, int j, int xdrop, int kmer_len, bool useHOPC, bool iRev, bool jRev, int cap) {

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

    bool reverse = (twin == seedV) || (useHOPC && ( iRev != jRev ));

    // use caps
    // if ( cap ) {
    //   if (!matchCaps(cap, reverse, row, col, i, j, useHOPC, kmer_len)) {
    //     longestExtensionScore.score = 0;
    //     longestExtensionScore.seed = seed;
    //     longestExtensionScore.strand = '0';
    //     return longestExtensionScore;
    //   }
    // }

    if ( reverse )
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
