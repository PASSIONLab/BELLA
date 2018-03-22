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

typedef Seed<Simple> TSeed;

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
pair<int,TSeed> seqanAlOne(std::string & row, std::string & col, int rlen, int i, int j, int dropFactor) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    int longestExtensionTemp;
    std::pair<int,TSeed> longestExtensionScore;

    seqH = row;
    seqV = col;

    TSeed seed1(i, j, i+KMER_LENGTH, j+KMER_LENGTH);
    seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
    seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        Dna5StringReverseComplement twinRead(seqH);
        i = rlen-i-KMER_LENGTH;
        TSeed seed1(i, j, i+KMER_LENGTH, j+KMER_LENGTH);

        /* Perform match extension */
        longestExtensionTemp = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());

    } else longestExtensionTemp = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());

    longestExtensionScore = make_pair(longestExtensionTemp, seed1);
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
pair<int,Seed<Simple>> seqanAlGen(std::string & row, std::string & col, int rlen, int i, int j, int l, int m, int dropFactor) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    std::pair<int,int> longestExtensionTemp;
    std::pair<int,Seed<Simple>> longestExtensionScore;

    seqH = row;
    seqV = col;

    TSeed seed1(i, j, i+KMER_LENGTH, j+KMER_LENGTH);
    TSeed seed2(l, m, l+KMER_LENGTH, m+KMER_LENGTH);
    seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
    seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        Dna5StringReverseComplement twinRead(seqH);

        i = rlen-i-KMER_LENGTH;
        l = rlen-l-KMER_LENGTH;

        TSeed seed1(i, j, i+KMER_LENGTH, j+KMER_LENGTH);
        TSeed seed2(l, m, l+KMER_LENGTH, m+KMER_LENGTH);

        /* Perform match extension */
        longestExtensionTemp.first = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        longestExtensionTemp.second = extendSeed(seed2, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());

    } else
    {
        longestExtensionTemp.first = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        longestExtensionTemp.second = extendSeed(seed2, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
    }

    if(longestExtensionTemp.first > longestExtensionTemp.second)
        longestExtensionScore = make_pair(longestExtensionTemp.first, seed1);
    else longestExtensionScore = make_pair(longestExtensionTemp.second, seed2);

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
 */
pair<int,TSeed> seqanAlOneAllKmer(std::string & row, std::string & col, int rlen, std::vector<std::pair<int,int>> vpos, int dropFactor) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    int longestExtensionTemp;
    std::pair<int,TSeed> longestExtensionScore;

    seqH = row;
    seqV = col;

    std::vector<std::pair<int,int>>::iterator it;
    it = vpos.begin();

    TSeed seed(it->first, it->second, it->first+KMER_LENGTH, it->second+KMER_LENGTH);
    seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
    seedV = infix(seqV, beginPositionV(seed), endPositionV(seed));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        Dna5StringReverseComplement twinRead(seqH);
        it->first = rlen-it->first-KMER_LENGTH;
        TSeed seed(it->first, it->second, it->first+KMER_LENGTH, it->second+KMER_LENGTH);

        /* Perform match extension */
        longestExtensionTemp = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());

    } else longestExtensionTemp = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());

    longestExtensionScore = make_pair(longestExtensionTemp, seed);
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
 */
pair<int,TSeed> seqanAlGenAllKmer(std::string & row, std::string & col, int rlen, std::vector<std::pair<int,int>> vpos, int dropFactor) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    int tempScore;
    std::pair<int,TSeed> longestExtensionScore;
    std::vector<std::pair<int,int>>::iterator fit;

    seqH = row;
    seqV = col;
    longestExtensionScore.first = 0;
    fit = vpos.begin();

    TSeed seed(fit->first, fit->second, fit->first+KMER_LENGTH, fit->second+KMER_LENGTH);
    seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
    seedV = infix(seqV, beginPositionV(seed), endPositionV(seed));
    
    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {   /* Reverse seq */
        Dna5StringReverseComplement twinRead(seqH);
        for(std::vector<std::pair<int,int>>::iterator it=vpos.begin(); it!=vpos.end(); ++it)
        {   /* Update position on reversed seq */
            it->first = rlen-it->first-KMER_LENGTH;
            /* Seed creation */
            TSeed seed(it->first, it->second, it->first+KMER_LENGTH, it->second+KMER_LENGTH);
            /* Perform match extension */
            tempScore = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
            /* Keep the best score */
            if(tempScore > longestExtensionScore.first)
                longestExtensionScore = make_pair(tempScore, seed);
        }
    } 
    else
    {
        for(std::vector<std::pair<int,int>>::iterator it=vpos.begin(); it!=vpos.end(); ++it)
        {   /* Seed creation */
            TSeed seed(it->first, it->second, it->first+KMER_LENGTH, it->second+KMER_LENGTH);
            /* Perform match extension */
            tempScore = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
            /* Keep the best score */
            if(tempScore > longestExtensionScore.first)
                longestExtensionScore = make_pair(tempScore, seed); 
        }  
    }
    
    return longestExtensionScore;
}
