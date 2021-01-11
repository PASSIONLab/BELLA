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

#include "common/common.h"

//	GG: check orientation
bool
checkstrand(const std::string& read1, const std::string& read2, const int begpH,
	const int begpV, const int kmerSize) {

	std::string seedH = read1.substr(begpH, kmerSize);
	std::string seedV = read2.substr(begpV, kmerSize);

	if(seedH != seedV) return false;
	else return true;
}

//	GG: check strand and compute overlap length
int
overlapop(const std::string& read1, const std::string& read2, unsigned short int begpH, 
	unsigned short int begpV, const unsigned short int kmerSize) {

	int read1len = read1.length();
	int read2len = read2.length();

	// GG: checking strand
	bool oriented = checkstrand(read1, read2, begpH, begpV, kmerSize);

	if(!oriented)
	{
		begpH = read1.length() - begpH - kmerSize;
	}

	// GG: computing overlap length
	unsigned short int endpH = begpH + kmerSize;
	unsigned short int endpV = begpV + kmerSize;

	int margin1 = std::min(begpH, begpV);
	int margin2 = std::min(read1len - endpH, read2len - endpV);
	int overlap = margin1 + margin2 + kmerSize;

	return overlap;
}

//	GG: multiply operation
void
multiop(spmatPtr_& value, const std::string& read1, const std::string& read2, 
	unsigned short int begpH, unsigned short int begpV, const int kmerSize) {

	value->count = 1;
	vector<pair<unsigned short int, unsigned short int>> vec{std::make_pair(begpH, begpV)};	// GG: to facilitate bin division
	value->pos.push_back(vec);
	value->support.push_back(1);	// initial k-mer has support 1

	// GG: check strand and compute overlap length
	int overlap = overlapop(read1, read2, begpH, begpV, kmerSize);
	value->overlap.push_back(overlap);
}

template <typename T,typename U>
std::pair<T,U> distance(const std::pair<T,U>& l, const std::pair<T,U>& r) {
	return {std::abs(l.first - r.first), std::abs(l.second - r.second)};
}

template <typename T,typename U>
bool operator>(const std::pair<T,U>& l, const T& c) {
	if(l.first > c && l.second > c) return true;
	else return false;
}

//	GG: binning operation
void
chainop(spmatPtr_& m1, spmatPtr_& m2, BELLApars& b_pars, 
	const std::string& readname1, const std::string& readname2)
{
	// number of common k-mer
	m1->count = m1->count + m2->count;
	vector<unsigned short int> tobeinserted;
	vector<vector<pair<unsigned short int, unsigned short int>>> kmertobeinserted(m1->pos.size());

	for(int i = 0; i < m2->pos.size(); ++i)	
	{
		bool orphan = true;
		for(int j = 0; j < m1->pos.size(); ++j)
		{
			if(std::abs(m2->overlap[i] - m1->overlap[j]) < b_pars.binSize)
			{
				for(auto kmer1:m1->pos[j])
				{
					for(auto kmer2:m2->pos[i])
					{
						//	GG: kmer need to be not overlapping and at least <kmerRift> distant from each other (kmerRift = kmerSize deafult)
						if(distance(kmer1, kmer2) > b_pars.kmerSize)
						{
							kmertobeinserted[j].push_back(kmer2);
						}
					}
				}
				orphan = false;	// we can be within b length of multiple overlap estimations, so we can't break
			}

			if(orphan)
			{
				tobeinserted.push_back(i);	// we don't want to immediately insert to m1 and increase computational complexity
			}
		}
	}

	for (int j = 0; j < kmertobeinserted.size(); j++)
	{
		m1->support[j] += kmertobeinserted[j].size();
		m1->count	   += kmertobeinserted[j].size();
		m1->pos[j].insert(m1->pos[j].end(), kmertobeinserted[j].begin(), kmertobeinserted[j].end());
	}

	for (auto i:tobeinserted)
	{
		m1->pos.push_back(m2->pos[i]);
		m1->overlap.push_back(m2->overlap[i]);
		m1->support.push_back(m2->support[i]);
	}
}

#endif
