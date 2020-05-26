#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/seeds.h>
#include "common.h"
#ifndef __NVCC__
	#include "../xavier/xavier.h"
#endif
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

#ifdef __NVCC__
#include "../loganGPU/logan.cuh"
#endif

#define BATCH_SIZE 30000

using namespace seqan;
using namespace std;

char
complementbase(char n) {
	switch(n)
	{
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	}
	assert(false);
	return ' ';
}

std::string
reversecomplement(const std::string& seq) {

	std::string cpyseq = seq;
	std::reverse(cpyseq.begin(), cpyseq.end());

	std::transform(
		std::begin(cpyseq),
		std::end  (cpyseq),
		std::begin(cpyseq),
	complementbase);

	return cpyseq;
}

double slope(double error)
{
	double p_mat = pow(1-error,2);  // match
	double p_mis = 1-p_mat;         // mismatch/gap
	double alpha = 1;               // match penalty
	double beta  = 1;                // mismatch/gap penalty

	return alpha*p_mat - beta*p_mis;
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

    bool reverse = (twin == seedV) || (useHOPC && ( iRev != jRev ));

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

#ifndef __NVCC__
/**
 * @brief alignLogan does the seed-and-extend alignment
 * @param row
 * @param col
 * @param rowLen is the length of the row sequence
 * @param i is the starting position of the k-mer on the first read
 * @param j is the starting position of the k-mer on the second read
 * @param xDrop
 * @return alignment score and extended seed
 */
xavierResult xavierAlign(const std::string& row, const std::string& col, int rowLen, int i, int j, int xDrop, int kmerSize)
{
	// result.first = best score, result.second = exit score when (if) x-drop termination is satified
	std::pair<int, int> tmp;
	xavierResult result;

	// penalties (LOGAN currently supports only linear gap penalty and penalty within +/- 3)
	short match    =  1;
	short mismatch = -1;
	short gap 	   = -1;

	// initialize scoring scheme
	ScoringSchemeX scoringScheme(match, mismatch, gap);	// enalties (LOGAN currently supports only linear gap penalty and penalty within +/- 3)

	SeedX seed(i, j, kmerSize);

	std::string seedH = row.substr(getBeginPositionH(seed), kmerSize);
	std::string seedV = col.substr(getBeginPositionV(seed), kmerSize);

	std::string rep = reversecomplement(seedH);
	std::string cpyrow(row);

	if(rep == seedV)
	{
		std::reverse(std::begin(cpyrow), std::end(cpyrow));
		std::transform(std::begin(cpyrow), std::end(cpyrow), std::begin(cpyrow), complementbase);

		setBeginPositionH(seed, rowLen - i - kmerSize);
		setEndPositionH(seed, rowLen - i);

		// perform match extension reverse string
 		tmp = XavierXDrop(seed, XAVIER_EXTEND_BOTH, cpyrow, col, scoringScheme, xDrop);
		result.strand = "c";
	}
	else
	{
		// perform match extension forward string
	 	tmp = XavierXDrop(seed, XAVIER_EXTEND_BOTH, row, col, scoringScheme, xDrop);
		result.strand = "n";
	}

	result.score = tmp.first; 	// best score

	setBeginPositionH(result.seed, getBeginPositionH(seed));	// updated extension
	setBeginPositionV(result.seed, getBeginPositionV(seed));	// updated extension

	setEndPositionH(result.seed, getEndPositionH(seed));		// updated extension
	setEndPositionV(result.seed, getEndPositionV(seed));		// updated extension

	return result;
}

#else

// ======================================= //
// 				GPU Functions			   //
// ======================================= //

void alignLogan(vector<string>&	target, vector<string>&	query, vector<SeedL>& seeds, 
	const BELLApars& b_pars, vector<loganResult>& longestExtensionScore)
{

	ScoringSchemeL sscheme(1, -1, -1, -1);
	std::vector<ScoringSchemeL> scoring;
	scoring.push_back(sscheme);

	int AlignmentsToBePerformed = seeds.size();
	printLog(AlignmentsToBePerformed);
	//int* res = (int*)malloc(BATCH_SIZE*sizeof(int));
	int numAlignmentsLocal = BATCH_SIZE * b_pars.numGPU; 
	cout <<"///////////////////////////////////////////////" <<b_pars.numGPU << endl;
	

	//	Divide the alignment in batches of 100K alignments
	for(int i = 0; i < AlignmentsToBePerformed; i += BATCH_SIZE * b_pars.numGPU)
	{
		if(AlignmentsToBePerformed < (i + BATCH_SIZE * b_pars.numGPU))
			numAlignmentsLocal = AlignmentsToBePerformed % (BATCH_SIZE * b_pars.numGPU);

		int* res = (int*)malloc(numAlignmentsLocal * sizeof(int));	

		std::vector<string>::const_iterator first_t = target.begin() + i;
		std::vector<string>::const_iterator last_t  = target.begin() + i + numAlignmentsLocal;
		std::vector<string> target_b(first_t, last_t);

		std::vector<string>::const_iterator first_q = query.begin() + i;
		std::vector<string>::const_iterator last_q  = query.begin() + i + numAlignmentsLocal;
		std::vector<string> query_b(first_q, last_q);

		std::vector<SeedL>::const_iterator first_s = seeds.begin() + i;
		std::vector<SeedL>::const_iterator last_s  = seeds.begin() + i + numAlignmentsLocal;
		std::vector<SeedL> seeds_b(first_s, last_s);

		extendSeedL(seeds_b, EXTEND_BOTHL, target_b, query_b, scoring, b_pars.xDrop, b_pars.kmerSize, res, numAlignmentsLocal, b_pars.numGPU);

		for(int j=0; j<numAlignmentsLocal; j++)
		{
			longestExtensionScore[j+i].score = res[j];
			longestExtensionScore[j+i].seed = seeds_b[j];
		}

		free(res);
	}
}
#endif // __NVCC__
#endif // _ALIGNMENT_H_
