#ifndef BELLA_CHAIN_H_
#define BELLA_CHAIN_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <array>
#include <tuple>
#include <queue>
#include <memory>
#include <stack>
#include <functional>
#include <cstring>
#include <string.h>
#include <math.h>
#include <cassert>
#include <ios>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <map>
#include <unordered_map>
#include <omp.h>

#include "mtspgemm2017/common.h"

// GG: check orientation
bool
checkstrand(const std::string& read1, const std::string& read2, const int begpH,
	const int begpV, const int kmerSize) {

	std::string seedH = read1.substr(begpH, kmerSize);
	std::string seedV = read2.substr(begpV, kmerSize);

	if(seedH != seedV) return false;
	else return true;
}

// GG: check strand and compute overlap length
int
overlapop(const std::string& read1, const std::string& read2, int begpH, 
	int begpV, const int kmerSize) {

	int read1len = read1.length();
	int read2len = read2.length();

	// GG: checking strand
	bool oriented = checkstrand(read1, read2, begpH, begpV, kmerSize);

	if(!oriented)
	{
		begpH = read1.length() - begpH - kmerSize;
	}

	// GG: computing overlap length
	int endpH = begpH + kmerSize;
	int endpV = begpV + kmerSize;

	int margin1 = std::min(begpH, begpV);
	int margin2 = std::min(read1len - endpH, read2len - endpV);
	int overlap = margin1 + margin2 + kmerSize;

	return overlap;
}

// GG: multiply operation
void
multiop(spmatPtr_& value, const std::string& read1, const std::string& read2, 
	int begpH, int begpV, const int kmerSize) {

	value->count = 1;
	value->pos.push_back(make_pair(begpH, begpV));
	value->support.push_back(1);	// initial k-mer has support 1

	// GG: check strand and compute overlap length
	int overlap = overlapop(read1, read2, begpH, begpV, kmerSize);
	value->overlap.push_back(overlap);
}

// GG: chaining operation
void
chainop(spmatPtr_& m1, spmatPtr_& m2, BELLApars& b_parameters, 
	const std::string& readname1, const std::string& readname2)
{

	// number of common k-mer
	m1->count = m1->count + m2->count;

	vector<int> tobeinserted;
	for(int i = 0; i < m2->pos.size(); ++i)	
	{
		bool orphan = true;
		for(int j = 0; j < m1->pos.size(); ++j)
		{
			// GG: TODO 500 as parameter
			if(std::abs(m2->overlap[i] - m1->overlap[j]) < b_parameters.bin) // B is the bin length
			{
				m1->support[j] += m2->support[j];
				orphan = false;
				// we can be within (B=500) length of multiple overlap estimations, so we can't break
			}
		}

		if(orphan)
		{
			tobeinserted.push_back(i);	// we don't want to immediately insert to m1 and increase computational complexity
		}
	}

	for (auto i:tobeinserted)
	{
		m1->pos.push_back(m2->pos[i]);
		m1->overlap.push_back(m2->overlap[i]);
		m1->support.push_back(m2->support[i]);
	}

#pragma omp critical
	{
		// GG: after testing correctness, this part can be removed
		// GG: we can then pass fewer parameters
		std::cout << "Between " << readname1 << " and " << readname2 << std::endl;
		std::copy(m1->overlap.begin(), m1->overlap.end(), std::ostream_iterator<int>(std::cout, " ")); std::cout << std::endl;
		std::copy(m1->support.begin(), m1->support.end(), std::ostream_iterator<int>(std::cout, " ")); std::cout << std::endl;
	}

}

#endif