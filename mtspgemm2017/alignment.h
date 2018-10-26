#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/seeds.h>
#include "global.h"
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

using namespace seqan;
using namespace std;

/**
 * @brief seqanAlOne does the seed-and-extend alignment
 * when ony one shared k-mer exists
 * @param row
 * @param col
 * @param rlen is the length of the row sequence
 * @param i is the starting position of the k-mer on the first read
 * @param j is the starting position of the k-mer on the second read
 * @param dropFactor
 * @return alignment score and extended seed
 */

seqAnResult seqanAlOne(std::string & row, std::string & col, int rlen, int i, int j, int dropFactor, int kmer_len) {

    Score<int, Simple> scoringScheme(1, -1, -1, -7);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    string strand;
    int longestExtensionTemp;
    seqAnResult longestExtensionScore;

    seqH = row;
    seqV = col;

    TSeed seed1(i, j, i+kmer_len, j+kmer_len);
    seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
    seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        strand = 'c';
        Dna5StringReverseComplement twinRead(seqH);
        i = rlen-i-kmer_len;
        
        setBeginPositionH(seed1, i);
        setBeginPositionV(seed1, j);
        setEndPositionH(seed1, i+kmer_len);
        setEndPositionV(seed1, j+kmer_len);

        /* Perform match extension */
        longestExtensionTemp = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());

    } else
    {
        strand = 'n';
        longestExtensionTemp = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
    } 

    longestExtensionScore.score = longestExtensionTemp;
    longestExtensionScore.seed = seed1;
    longestExtensionScore.strand = strand;
    return longestExtensionScore;
}
/**
 * @brief seqanAlGen does the seed-and-extend alignment
 * when two shared k-mers exist
 * @param row
 * @param col
 * @param rlen is the length of the row sequence
 * @param i is the starting position of the first k-mer on the first read
 * @param j is the starting position of the first k-mer on the second read
 * @param l is the starting position of the second k-mer on the first read
 * @param m is the starting position of the second k-mer on the second read
 * @param dropFactor
 * @return alignment score and extended seed
 */
seqAnResult seqanAlGen(std::string & row, std::string & col, int rlen, int i, int j, int l, int m, int dropFactor, int kmer_len) {

    Score<int, Simple> scoringScheme(1, -1, -1, -7);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    string strand;
    std::pair<int,int> longestExtensionTemp;
    seqAnResult longestExtensionScore;

    seqH = row;
    seqV = col;

    TSeed seed1(i, j, i+kmer_len, j+kmer_len);
    TSeed seed2(l, m, l+kmer_len, m+kmer_len);
    seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
    seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        strand = 'c';
        Dna5StringReverseComplement twinRead(seqH);

        i = rlen-i-kmer_len;
        l = rlen-l-kmer_len;

        setBeginPositionH(seed1, i);
        setBeginPositionV(seed1, j);
        setEndPositionH(seed1, i+kmer_len);
        setEndPositionV(seed1, j+kmer_len);

        setBeginPositionH(seed2, l);
        setBeginPositionV(seed2, m);
        setEndPositionH(seed2, l+kmer_len);
        setEndPositionV(seed2, m+kmer_len);

        /* Perform match extension */
        longestExtensionTemp.first = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        longestExtensionTemp.second = extendSeed(seed2, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());

    } else
    {
        strand = 'n';
        longestExtensionTemp.first = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        longestExtensionTemp.second = extendSeed(seed2, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
    }

    longestExtensionScore.strand = strand;

    if(longestExtensionTemp.first > longestExtensionTemp.second)
    {
        longestExtensionScore.score = longestExtensionTemp.first;
        longestExtensionScore.seed = seed1;
    }
    else
    {
        longestExtensionScore.score = longestExtensionTemp.second;
        longestExtensionScore.seed = seed2;
    } 

    return longestExtensionScore;
}
/**
 * @brief seqanAlOneAllKmer does the seed-and-extend alignment
 * when only one shared k-mer exists
 * @param row
 * @param col
 * @param rlen is the length of the row sequence
 * @param vpos is a vector containing the starting position of k-mer
 * @param dropFactor
 * @return alignment score and extended seed
 * DO NOT USE IT UNTIL THE SEED POSITION DEFINITION WILL BE CORRECTED AND UPDATED
 */
seqAnResult seqanAlOneAllKmer(std::string & row, std::string & col, int rlen, std::vector<std::pair<int,int>> vpos, int dropFactor, int kmer_len) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    string strand;
    int longestExtensionTemp;
    seqAnResult longestExtensionScore;

    seqH = row;
    seqV = col;

    std::vector<std::pair<int,int>>::iterator it;
    it = vpos.begin();

    TSeed seed(it->first, it->second, it->first+kmer_len, it->second+kmer_len);
    seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
    seedV = infix(seqV, beginPositionV(seed), endPositionV(seed));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        strand = 'c';
        Dna5StringReverseComplement twinRead(seqH);
        it->first = rlen-it->first-kmer_len;

        setBeginPositionH(seed, it->first);
        setBeginPositionV(seed, it->second);
        setEndPositionH(seed, it->first+kmer_len);
        setEndPositionV(seed, it->second+kmer_len);

        //TSeed seed(it->first, it->second, it->first+kmer_len, it->second+kmer_len); // In this way the seed is not globally updated

        /* Perform match extension */
        longestExtensionTemp = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());

    } else
    {
        strand = 'n';
        longestExtensionTemp = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
    } 

    longestExtensionScore.score = longestExtensionTemp;
    longestExtensionScore.seed = seed;
    longestExtensionScore.strand = strand;

    return longestExtensionScore;
}

/**
 * @brief seqanAlGenAllKmer does the seed-and-extend alignment
 * when shared k-mers > 1
 * @param row
 * @param col
 * @param rlen
 * @param vpos
 * @param dropFactor
 * @return alignment score and extended seed
 * DO NOT USE IT UNTIL THE SEED POSITION DEFINITION WILL BE CORRECTED AND UPDATED
 */
seqAnResult seqanAlGenAllKmer(std::string & row, std::string & col, int rlen, std::vector<std::pair<int,int>> vpos, int dropFactor, int kmer_len) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    string strand;
    int tempScore;
    seqAnResult longestExtensionScore;
    std::vector<std::pair<int,int>>::iterator fit;

    seqH = row;
    seqV = col;
    longestExtensionScore.score = 0;
    fit = vpos.begin();

    TSeed seed(fit->first, fit->second, fit->first+kmer_len, fit->second+kmer_len);
    seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
    seedV = infix(seqV, beginPositionV(seed), endPositionV(seed));
    
    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {   /* Reverse seq */
        strand = 'c';
        Dna5StringReverseComplement twinRead(seqH);
        for(std::vector<std::pair<int,int>>::iterator it=vpos.begin(); it!=vpos.end(); ++it)
        {   /* Update position on reversed seq */
            it->first = rlen-it->first-kmer_len;
            /* Seed unpdate */
            setBeginPositionH(seed, it->first);
            setBeginPositionV(seed, it->second);
            setEndPositionH(seed, it->first+kmer_len);
            setEndPositionV(seed, it->second+kmer_len);
            //TSeed seed(it->first, it->second, it->first+kmer_len, it->second+kmer_len); // In this way the seed is not globally updated
            /* Perform match extension */
            tempScore = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
            /* Keep the best score */
            if(tempScore > longestExtensionScore.score)
            {
                longestExtensionScore.score = tempScore;
                longestExtensionScore.seed = seed;
            }
        }
    } 
    else
    {
        strand = 'n';
        for(std::vector<std::pair<int,int>>::iterator it=vpos.begin(); it!=vpos.end(); ++it)
        {   /* Seed update */
            setBeginPositionH(seed, it->first);
            setBeginPositionV(seed, it->second);
            setEndPositionH(seed, it->first+kmer_len);
            setEndPositionV(seed, it->second+kmer_len);
            //TSeed seed(it->first, it->second, it->first+kmer_len, it->second+kmer_len);
            /* Perform match extension */
            tempScore = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
            /* Keep the best score */
            if(tempScore > longestExtensionScore.score)
            {
                longestExtensionScore.score = tempScore;
                longestExtensionScore.seed = seed; 
            }
        }  
    }
    longestExtensionScore.strand = strand;

    return longestExtensionScore;
}

#endif
