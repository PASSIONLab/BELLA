#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <map>
#include <sys/types.h> 
#include <sys/stat.h> 
#include <math.h>
#include <limits.h>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h> 
#include <sstream>
#include <set>
#include <memory>
#include <typeinfo>
#include <set>
#include <iomanip>
#include "def.h"

std::multiset<entry, classcom> cutShortOverlaps(std::multiset<entry, classcom>& s, int minOverlap) {

	std::multiset<entry, classcom> Sset;

	for(auto it = s.begin(); it != s.end(); it++)
		if(it->overlap >= minOverlap)
			Sset.insert(*it);

	return Sset;
}

void estimateOverlap(int begV, int endV, int lenV, int begH, int endH, int lenH, int& overlap)
{
	overlap = min(begV, begH) + min(lenV - endV, lenH - endH) + ((endV - begV) + (endH - begH)) / 2;
}

std::multiset<entry, classcom> readTruthOutput(std::ifstream& file, int minOverlap, bool isSimulated)
{
	int maxt = 1;
#pragma omp parallel
	{
		maxt = omp_get_num_threads();
	}

	int nOverlap = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
	file.seekg(0,std::ios_base::beg);

	std::vector<std::string> entries;
	entries.reserve(nOverlap);

	if(file)
		for (int i = 0; i < nOverlap; ++i) {
			std::string line;
			std::getline(file, line);
			entries.push_back(line);
		}
	file.close();

	std::map<std::string, std::vector<overlap>> intermediate;
	std::vector<std::map<std::string, std::vector<overlap>>> local(maxt);

#pragma omp parallel for
	for(int i = 0; i < nOverlap; i++) {

		int ithread = omp_get_thread_num();

		std::vector<std::string> v = split(entries[i], ' ');
		overlap ioverlap;

		if(isSimulated) {
			ioverlap.ref   =       v[0];
			ioverlap.read  = "@" + v[3];
			ioverlap.start =       stoi(v[1]);
			ioverlap.end   =       stoi(v[2]);
		} else {
			ioverlap.ref   =       v[0];
			ioverlap.read  = "@" + v[1];
			ioverlap.start =       stoi(v[2]);
			ioverlap.end   =       stoi(v[3]);
		}
		auto it = local[ithread].find(ioverlap.ref);

		if(it == local[ithread].end()) {
			std::vector<overlap> tmp;
			tmp.push_back(ioverlap);
			local[ithread].insert(std::map<std::string, std::vector<overlap>>::value_type(ioverlap.ref, tmp));
		}
		else {
			it->second.push_back(ioverlap);
			local[ithread][ioverlap.ref] = it->second;
		}
	}

	for(int i = 0; i < maxt; ++i)
	{
		intermediate = std::accumulate(local[i].begin(), local[i].end(), intermediate,
				[](std::map<std::string, std::vector<overlap>>& m, 
					const std::pair<std::string, std::vector<overlap>>& v)
				{
					m[v.first].insert(m[v.first].end(), v.second.begin(), v.second.end());
					return m;
				});
	}

	vector<Interval<std::string>> intervals;
	vector<Interval<std::string>> queries;
	std::multiset<entry, classcom> Gset;

	for(auto key = intermediate.begin(); key != intermediate.end(); ++key) {
		for(auto it = key->second.begin(); it != key->second.end(); ++it)
		{
			intervals.push_back(Interval<std::string>(it->start, it->end, it->read));
			queries.push_back(Interval<std::string>(it->start, it->end, it->read));
		}

		IntervalTree<string> tree;
		tree = IntervalTree<std::string>(intervals);

		for (auto q = queries.begin(); q != queries.end(); ++q) {
			vector<Interval<std::string>> results;
			tree.findOverlapping(q->start, q->stop, q->value, results, minOverlap, Gset);
		}

		intervals.clear();
		queries.clear();
	}

#ifdef DEBUG
	std::cout << std::endl;
	std::cout << Gset.size() << " overlaps in the ground truth longer than " << minOverlap << " bp"<< std::endl;
	std::cout << std::endl;
#endif

	return Gset;
};

std::multiset<entry, classcom> readBellaOutput(std::ifstream& file)
{
	int maxt = 1;
#pragma omp parallel
	{
		maxt = omp_get_num_threads();
	}

	int nOverlap = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
	file.seekg(0,std::ios_base::beg);
	std::vector<std::string> entries;

	if(file)
		for (int i = 0; i < nOverlap; ++i) {
			std::string line;
			std::getline(file, line);
			entries.push_back(line);
		}
	file.close();

	std::multiset<entry, classcom> result;
	std::vector<std::multiset<entry, classcom>> local(maxt);

#pragma omp parallel for
	for(int i = 0; i < nOverlap; i++) {
		std::stringstream linestream(entries[i]);
		int ithread = omp_get_thread_num();

		std::vector<std::string> v = split(entries[i], '\t');
		entry ientry;

		//	std::cout << "What's up, dude?" << std::endl;
		ientry.a = v[0];
		ientry.b = v[1];
		ientry.overlap = stoi(v[4]);

		local[ithread].insert(ientry);
	}

	for(int i = 0; i < maxt; ++i)
		result.insert(local[i].begin(), local[i].end());
#ifdef DEBUG
	std::cout << "BELLA identified " << 2*result.size() << " overlaps." << std::endl;
#endif
	return result;
};

std::multiset<entry, classcom> readMinimapOutput(std::ifstream& file)
{
	int maxt = 1;
#pragma omp parallel
	{
		maxt = omp_get_num_threads();
	}

	int nOverlap = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
	file.seekg(0,std::ios_base::beg);
	std::vector<std::string> entries;

	if(file)
		for (int i = 0; i < nOverlap; ++i) {
			std::string line;
			std::getline(file, line);
			entries.push_back(line);
		}
	file.close();

	std::multiset<entry, classcom> result;
	std::vector<std::multiset<entry, classcom>> local(maxt);

#pragma omp parallel for
	for(int i = 0; i < nOverlap; i++) {
		std::stringstream linestream(entries[i]);
		int ithread = omp_get_thread_num();

		std::vector<std::string> v = split(entries[i], '\t');
		entry ientry;

	//	std::cout << "What's up, dude?" << std::endl;
		ientry.a = "@" + v[0];
		ientry.b = "@" + v[5];

	//	paf format
		int lenV = std::stoi(v[1]);
		int begV = std::stoi(v[2]);
		int endV = std::stoi(v[3]);
		int lenH = std::stoi(v[6]);
		int begH = std::stoi(v[7]);
		int endH = std::stoi(v[8]);

	//	coordinate are reported on the original strand in paf format
		if(v[4] == "-") {
			int tmp = begH;
			begH = lenH - endH;
			endH = lenH - tmp;
		}

	//	(int begV, int endV, int lenV, int begH, int endH, int lenH, int overlap)
		estimateOverlap(begV, endV, lenV, 
				begH, endH, lenH, ientry.overlap);

		local[ithread].insert(ientry);
	}

	for(int i = 0; i < maxt; ++i)
		result.insert(local[i].begin(), local[i].end());
#ifdef DEBUG
	std::cout << "MINIMAP2 identified " << 2*result.size() << " overlaps." << std::endl;
#endif
	return result;
};

//std::set<entry, classcom> readMecatOutput(std::ifstream& output)
//{
//
//};
//
//std::set<entry, classcom> readMhapOutput(std::ifstream& output)
//{
//
//};
//
//std::set<entry, classcom> readBlasrOutput(std::ifstream& output)
//{
//
//};
//
//std::set<entry, classcom> readDalignerOutput(std::ifstream& output)
//{
//
//};

void evaluate(std::multiset<entry, classcom>& Sset, const std::multiset<entry, classcom>& Gset, int minOverlap, bool toDuplicate)
{
	std::multiset<entry, classcom> Tset;
	Sset = cutShortOverlaps(Sset, minOverlap);

	std::set_intersection(
		Gset.begin(), Gset.end(),
		Sset.begin(), Sset.end(),
		std::inserter(Tset, Tset.end())
	);

	float RC;
	if(toDuplicate)
	{
	#ifdef DEBUG
		std::cout << "\t* " << 2*Sset.size() << " overlaps longer than " << minOverlap << " bp" << std::endl;
		std::cout << "\t* " << 2*Tset.size() << " true positives" << std::endl;
	#endif
		RC   = ((float)(2 * Tset.size()) / (float)Gset.size()) * 100;
	}
	else
	{
	#ifdef DEBUG
		std::cout << "\t* " << Sset.size() << " overlaps longer than " << minOverlap << " bp" << std::endl;
		std::cout << "\t* " << Tset.size() << " true positives" << std::endl;
	#endif
		RC   = ((float)Tset.size()       / (float)Gset.size()) * 100;
	}

	float PR = ((float)Tset.size()       / (float)Sset.size()) * 100;
	float F1 = (2 * RC * PR) / (RC + PR);

	std::cout << std::setprecision(2) 	<< std::fixed;
	std::cout << "Recall\t\t" 	<< RC << "%"<< std::endl;
	std::cout << "Precision\t" 	<< PR << "%"<< std::endl;
	std::cout << "F1\t\t" 		<< F1 << "%"<< std::endl;
	std::cout << std::endl;
}
