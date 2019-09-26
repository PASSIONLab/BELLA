#ifndef BELLA_KMERCOUNT_H_
#define BELLA_KMERCOUNT_H_

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
#include <numeric>
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
#include <fstream>
#include <typeinfo>

#ifdef __NVCC__
#include "loganGPU/logan.cuh"
#endif

#include "libcuckoo/cuckoohash_map.hh"
#include "libbloom/bloom64.h"

#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include "kmercode/Buffer.h"
#include "kmercode/common.h"
#include "kmercode/fq_reader.h"
#include "kmercode/ParallelFASTQ.h"
#include "kmercode/bound.hpp"
#include "kmercode/hyperloglog.hpp"
#include "mtspgemm2017/common.h"

using namespace std;

#define ASCIIBASE 33 // Pacbio quality score ASCII BASE

#ifndef PRINT
#define PRINT
#endif

//	GG: when couting k-mers the values can be 16bit (k-mer occurrence) instead of 32 (k-mer ids in the final dictionary)
typedef cuckoohash_map<Kmer, unsigned int>       dictionary_t_32bit;	// <k-mer && reverse-complement, #kmers>
typedef cuckoohash_map<Kmer, unsigned short int> dictionary_t_16bit;	// <k-mer && reverse-complement, kmer_id>

struct filedata {

	char filename[MAX_FILE_PATH];
	size_t filesize;
};

/**
 * @brief GetFiles
 * @param filename
 * @return
 */
vector<filedata>  GetFiles(char *filename) {
	int64_t totalsize = 0;
	int numfiles = 0;
	std::vector<filedata> filesview;
	
	filedata fdata;
	ifstream allfiles(filename);
	if(!allfiles.is_open()) {
		std::string someString(filename);
		std::string ErrorMessage = "Could not open " + someString;
		printLog(ErrorMessage);
		exit(1);
	}
	allfiles.getline(fdata.filename, MAX_FILE_PATH);
	while(!allfiles.eof())
	{
		struct stat st;
		stat(fdata.filename, &st);
		fdata.filesize = st.st_size;
		
		filesview.push_back(fdata);
		std::string InputFile = filesview.back().filename;
		std::string str1 = std::to_string((float)filesview.back().filesize / (1024*1024));
		std::string str2 = " MB";
		std::string InputSize = str1 + str2;

		printLog(InputFile);
		printLog(InputSize);
		allfiles.getline(fdata.filename, MAX_FILE_PATH);
		totalsize += fdata.filesize;
		numfiles++;
	}
	return filesview;
}

/**
 * @brief SimpleCount
 * @param allfiles
 * @param countsreliable_denovo
 * @param lower
 * @param upper
 * @param b_pars.kmerSize
 * @param upperlimit
 */
void SimpleCount(vector<filedata> & allfiles, dictionary_t_32bit& countsreliable_denovo, int& lower, int& upper,
	int coverage, size_t upperlimit, BELLApars & b_pars)
{
	vector < vector<double> > allquals(MAXTHREADS);

	double denovocount = omp_get_wtime();
	size_t totreads = 0;

	dictionary_t_16bit countsdenovo;
	auto updatefn = [](unsigned short int &count) { if (count < std::numeric_limits<unsigned short int>::max()) ++count; };

	for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) 
	{
		#pragma omp parallel
		{
			ParallelFASTQ *pfq = new ParallelFASTQ();
			pfq->open(itr->filename, false, itr->filesize);

			if(MYTHREAD == 0)
    			{
    		    		const char* ReadingFASTQ = itr->filename;
    		    		printLog(ReadingFASTQ);
    			}

			vector<string> seqs;
			vector<string> quals;
			vector<string> nametags;
			size_t tlreads = 0; // thread local reads

			size_t fillstatus = 1;
			while(fillstatus) 
			{ 
				fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
				size_t nreads = seqs.size();

				for(int i=0; i<nreads; i++) 
				{
					// remember that the last valid position is length()-1
					int len = seqs[i].length();
					double rerror = 0.0;

					for(int j = 0; j<= len - b_pars.kmerSize; j++)  
					{
						std::string kmerstrfromfastq = seqs[i].substr(j, b_pars.kmerSize);
						Kmer mykmer(kmerstrfromfastq.c_str(), kmerstrfromfastq.length());
						Kmer lexsmall = mykmer.rep();

                        			countsdenovo.upsert(lexsmall, updatefn, 1);
						
						if(b_pars.skipEstimate == false) 
						{
								// accuracy
								int bqual = (int)quals[i][j] - ASCIIBASE;
								double berror = pow(10,-(double)bqual/10);
								rerror += berror;
						}
					}

					if(b_pars.skipEstimate == false) 
					{
						// remaining k qual position accuracy
						for(int j = len - b_pars.kmerSize + 1; j < len; j++)
						{
							int bqual = (int)quals[i][j] - ASCIIBASE;
							double berror = pow(10,-(double)bqual/10);
							rerror += berror;
						}
						rerror = rerror / len;
						allquals[MYTHREAD].push_back(rerror);
					}
				} // for(int i=0; i<nreads; i++)
				tlreads += nreads;
			} //while(fillstatus) 
			delete pfq;

			#pragma omp critical
			{
				totreads += tlreads;
			}
		}
	}

	// Error estimation
	double& errorRate = b_pars.errorRate;
	if(b_pars.skipEstimate == false)
	{
		errorRate = 0.0; // reset to 0 here, otherwise it cointains default or user-defined values
		#pragma omp for reduction(+:errorRate)
		for (int i = 0; i < MAXTHREADS; i++) 
			{
				double temp = std::accumulate(allquals[i].begin(), allquals[i].end(), 0.0);
				errorRate += temp / (double)allquals[i].size();
			}
		b_pars.errorRate = errorRate / (double)MAXTHREADS;
	}

	double load2kmers = omp_get_wtime(); 
	std::string kmerCountingTime = std::to_string(load2kmers - denovocount) + " seconds";
	printLog(kmerCountingTime);


	
	// Reliable bounds computation using estimated error rate from phred quality score
	lower = computeLower(coverage, b_pars.errorRate, b_pars.kmerSize, b_pars.minProbability);
	upper = computeUpper(coverage, b_pars.errorRate, b_pars.kmerSize, b_pars.minProbability);

	// Reliable k-mer filter on countsdenovo
	unsigned int kmer_id_denovo = 0;
	auto lt = countsdenovo.lock_table(); // our counting
	for (const auto &it : lt) 
		if (it.second >= lower && it.second <= upper)
		{
			countsreliable_denovo.insert(it.first, kmer_id_denovo);
			++kmer_id_denovo;
		}
	lt.unlock(); // unlock the table

	// Print some information about the table
	if (countsreliable_denovo.size() == 0)
	{
		std::string ErrorMessage = "BELLA terminated: 0 entries within reliable range. You may want to reduce the k-mer lenght.";
		printLog(ErrorMessage);
		exit(1);
	} 
	else 
	{
		int numReliableKmers = countsreliable_denovo.size();
		printLog(numReliableKmers);
	}

	countsdenovo.clear(); // free
}

/**
 * @brief DeNovoCount
 * @param allfiles
 * @param countsreliable_denovo
 * @param lower
 * @param upper
 * @param b_pars.kmerSize
 * @param upperlimit
 */
void DeNovoCount(vector<filedata> & allfiles, dictionary_t_32bit& countsreliable_denovo, int& lower, int& upper,
	int coverage, size_t upperlimit, BELLApars & b_pars)
{
	vector < vector<Kmer> >   allkmers(MAXTHREADS);
	vector < vector<double> > allquals(MAXTHREADS);
	vector < HyperLogLog > hlls(MAXTHREADS, HyperLogLog(12));   // std::vector fill constructor

	double denovocount = omp_get_wtime();
	double CardinalityEstimate;
	size_t totreads = 0;

	for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) 
	{
		#pragma omp parallel
		{
			ParallelFASTQ *pfq = new ParallelFASTQ();
			pfq->open(itr->filename, false, itr->filesize);

			if(MYTHREAD == 0)
    			{
    		    		const char* ReadingFASTQ = itr->filename;
    		    		printLog(ReadingFASTQ);
    			}

			vector<string> seqs;
			vector<string> quals;
			vector<string> nametags;
			size_t tlreads = 0; // thread local reads

			size_t fillstatus = 1;
			while(fillstatus) 
			{ 
				fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
				size_t nreads = seqs.size();

				for(int i=0; i<nreads; i++) 
				{
					// remember that the last valid position is length()-1
					int len = seqs[i].length();
					double rerror = 0.0;

					for(int j = 0; j<= len - b_pars.kmerSize; j++)  
					{
						std::string kmerstrfromfastq = seqs[i].substr(j, b_pars.kmerSize);
						Kmer mykmer(kmerstrfromfastq.c_str(), kmerstrfromfastq.length());
						Kmer lexsmall = mykmer.rep();

						allkmers[MYTHREAD].push_back(lexsmall);
						hlls[MYTHREAD].add((const char*) lexsmall.getBytes(), lexsmall.getNumBytes());

						if(b_pars.skipEstimate == false) 
						{
								// accuracy
								int bqual = (int)quals[i][j] - ASCIIBASE;
								double berror = pow(10,-(double)bqual/10);
								rerror += berror;
						}
					}

					if(b_pars.skipEstimate == false) 
					{
						// remaining k qual position accuracy
						for(int j = len - b_pars.kmerSize + 1; j < len; j++)
						{
							int bqual = (int)quals[i][j] - ASCIIBASE;
							double berror = pow(10,-(double)bqual/10);
							rerror += berror;
						}
						rerror = rerror / len;
						allquals[MYTHREAD].push_back(rerror);
					}
				} // for(int i=0; i<nreads; i++)
				tlreads += nreads;
			} //while(fillstatus) 
			delete pfq;

			#pragma omp critical
			{
				totreads += tlreads;
			}
		}
	}

	// Error estimation
	double& errorRate = b_pars.errorRate;
	if(b_pars.skipEstimate == false)
	{
		errorRate = 0.0; // reset to 0 here, otherwise it cointains default or user-defined values
		#pragma omp for reduction(+:errorRate)
		for (int i = 0; i < MAXTHREADS; i++) 
			{
				double temp = std::accumulate(allquals[i].begin(), allquals[i].end(), 0.0);
				errorRate += temp / (double)allquals[i].size();
			}
		b_pars.errorRate = errorRate / (double)MAXTHREADS;
	}

	// HLL reduction (serial for now) to avoid double iteration
	for (int i = 1; i < MAXTHREADS; i++) 
	{
		std::transform(hlls[0].M.begin(), hlls[0].M.end(), hlls[i].M.begin(), hlls[0].M.begin(), [](uint8_t c1, uint8_t c2) -> uint8_t{ return std::max(c1, c2); });
	}
	CardinalityEstimate = hlls[0].estimate();

	double load2kmers = omp_get_wtime(); 
	std::string kmerCountingTime = std::to_string(load2kmers - denovocount) + " seconds";
	printLog(kmerCountingTime);

	const double desired_probability_of_false_positive = 0.05;
	struct bloom * bm = (struct bloom*) malloc(sizeof(struct bloom));
	bloom_init64(bm, CardinalityEstimate * 1.1, desired_probability_of_false_positive);

	double TableSize = ((double)bm->bits)/8/1024/1024;
	int numHashFunctions = bm->hashes;

	printLog(CardinalityEstimate);
	printLog(TableSize);
	printLog(numHashFunctions);

	dictionary_t_16bit countsdenovo;

	#pragma omp parallel
	{
		for(auto v:allkmers[MYTHREAD])
		{
			bool inBloom = (bool) bloom_check_add(bm, v.getBytes(), v.getNumBytes(),1);
			if(inBloom) countsdenovo.insert(v, 0);
		}
	}

	double firstpass = omp_get_wtime();
	std::string FirstKmerPassTime = std::to_string(firstpass - load2kmers) + " seconds";
	printLog(FirstKmerPassTime);

	free(bm); // release bloom filter memory

	// in this pass, only use entries that already are in the hash table
	auto updatecount = [](unsigned short int &num) { ++num; };
	#pragma omp parallel
	{
		for(auto v:allkmers[MYTHREAD])
		{
			// does nothing if the entry doesn't exist in the table
			countsdenovo.update_fn(v, updatecount);
		}
	}

	std::string SecondKmerPassTime = std::to_string(omp_get_wtime() - firstpass) + " seconds";
	printLog(SecondKmerPassTime);

	// Reliable bounds computation using estimated error rate from phred quality score
	lower = computeLower(coverage, b_pars.errorRate, b_pars.kmerSize, b_pars.minProbability);
	upper = computeUpper(coverage, b_pars.errorRate, b_pars.kmerSize, b_pars.minProbability);

	// Reliable k-mer filter on countsdenovo
	unsigned int kmer_id_denovo = 0;
	auto lt = countsdenovo.lock_table(); // our counting
	for (const auto &it : lt) 
		if (it.second >= lower && it.second <= upper)
		{
			countsreliable_denovo.insert(it.first, kmer_id_denovo);
			++kmer_id_denovo;
		}
	lt.unlock(); // unlock the table

	// Print some information about the table
	if (countsreliable_denovo.size() == 0)
	{
		std::string ErrorMessage = "BELLA terminated: 0 entries within reliable range. You may want to reduce the k-mer lenght.";
		printLog(ErrorMessage);
		exit(1);
	} 
	else 
	{
		int numReliableKmers = countsreliable_denovo.size();
		printLog(numReliableKmers);
	}

	countsdenovo.clear(); // free
}

// Returns the new average after including x 
double getAvg(double prev_avg, int64_t x, int64_t n) 
{ 
	return (prev_avg * n + x) / (n + 1); 
} 

/**
 * @brief Split4Count
 * @param allfiles
 * @param countsreliable_denovo
 * @param lower
 * @param upper
 * @param b_pars.kmerSize
 * @param upperlimit
 */
void Split4Count(vector<filedata> & allfiles, dictionary_t_32bit& countsreliable_denovo, int& lower, int& upper, int coverage, size_t upperlimit, BELLApars & b_pars)
{
	size_t totreads = 0;
	size_t totbases = 0;
	
	double avesofar = 0.0;
	
	// Reliable k-mer filter on countsdenovo
	unsigned int kmer_id_denovo = 0;

	for(int splits=0; splits< 4; ++splits)	// splits
	{
		double denovocount = omp_get_wtime();
		
		vector < vector<Kmer> >   allkmers(MAXTHREADS);
		vector < HyperLogLog > hlls(MAXTHREADS, HyperLogLog(12));   // std::vector fill constructor
		
		for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) 
		{
			#pragma omp parallel
			{
				double tlave = 0.0;

				ParallelFASTQ *pfq = new ParallelFASTQ();
				pfq->open(itr->filename, false, itr->filesize);

				if(MYTHREAD == 0)
    				{
    		    			const char* ReadingFASTQ = itr->filename;
    		    			printLog(ReadingFASTQ);
    				}

				vector<string> seqs;
				vector<string> quals;
				vector<string> nametags;
				size_t tlreads = 0; // thread local reads
				size_t tlbases = 0; // thread local bases

				size_t fillstatus = 1;
				while(fillstatus) 
				{ 
					fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
					size_t nreads = seqs.size();

					for(int i=0; i<nreads; i++) 
					{
						// remember that the last valid position is length()-1
						int len = seqs[i].length();
						double rerror = 0.0;

						for(int j = 0; j<= len - b_pars.kmerSize; j++)  
						{
							std::string kmerstrfromfastq = seqs[i].substr(j, b_pars.kmerSize);
							Kmer mykmer(kmerstrfromfastq.c_str(), kmerstrfromfastq.length());
							Kmer lexsmall = mykmer.rep();
							if(lexsmall.hash() % 4 ==  splits)	// mod 4
							{
								allkmers[MYTHREAD].push_back(lexsmall);
								hlls[MYTHREAD].add((const char*) lexsmall.getBytes(), lexsmall.getNumBytes());

								if(b_pars.skipEstimate == false) 
								{
									// accuracy
									int bqual = (int)quals[i][j] - ASCIIBASE;
									double berror = pow(10,-(double)bqual/10);
        								tlave = getAvg(tlave, berror, tlbases++); 	
								}
							}
						}

						if(b_pars.skipEstimate == false) 
						{
							// remaining k qual position accuracy
							for(int j = len - b_pars.kmerSize + 1; j < len; j++)
							{
								int bqual = (int)quals[i][j] - ASCIIBASE;
								double berror = pow(10,-(double)bqual/10);
        							tlave = getAvg(tlave, berror, tlbases++); 									
							}
						}
					} // for(int i=0; i<nreads; i++)
					tlreads += nreads;
				} //while(fillstatus) 
				delete pfq;

				#pragma omp critical
				{
					totreads += tlreads;
					avesofar = (avesofar * totbases + tlave * tlbases) / (totbases + tlbases);				
					totbases += tlbases;
				}

			} // #pragma omp parallel

		} // for allfiles

		if(b_pars.skipEstimate == false)
		{
			b_pars.errorRate = avesofar;	// this 25% sample should be good enough, no?
			printLog(b_pars.errorRate);
		}
	
		// HLL reduction (serial for now) to avoid double iteration
		for (int i = 1; i < MAXTHREADS; i++) 
		{
			std::transform(hlls[0].M.begin(), hlls[0].M.end(), hlls[i].M.begin(), hlls[0].M.begin(), [](uint8_t c1, uint8_t c2) -> uint8_t{ return std::max(c1, c2); });
		}
		double CardinalityEstimate = hlls[0].estimate();

		double load2kmers = omp_get_wtime(); 
		std::string kmerCountingTime = std::to_string(load2kmers - denovocount) + " seconds";
		printLog(kmerCountingTime);
	
	
		const double desired_probability_of_false_positive = 0.05;
		struct bloom * bm = (struct bloom*) malloc(sizeof(struct bloom));
		bloom_init64(bm, CardinalityEstimate * 1.1, desired_probability_of_false_positive);

		double TableSize = ((double)bm->bits)/8/1024/1024;
		int numHashFunctions = bm->hashes;

		printLog(CardinalityEstimate);
		printLog(TableSize);
		printLog(numHashFunctions);

		dictionary_t_16bit countsdenovo;

		#pragma omp parallel
		{
			for(auto v:allkmers[MYTHREAD])
			{
				bool inBloom = (bool) bloom_check_add(bm, v.getBytes(), v.getNumBytes(),1);
				if(inBloom) countsdenovo.insert(v, 0);
			}
		}
		size_t tot_kmers = 0;
		for (int i=0; i<MAXTHREADS; i++)
			tot_kmers+= allkmers[i].size();

		printLog(splits);
		printLog(tot_kmers);	
	

		double firstpass = omp_get_wtime();
		std::string FirstKmerPassTime = std::to_string(firstpass - load2kmers) + " seconds";
		printLog(FirstKmerPassTime);

		free(bm); // release bloom filter memory

		// in this pass, only use entries that already are in the hash table
		auto updatecount = [](unsigned short int &num) { ++num; };
		#pragma omp parallel
		{
			for(auto v:allkmers[MYTHREAD])
			{
				// does nothing if the entry doesn't exist in the table
				countsdenovo.update_fn(v, updatecount);
			}
		}

		std::string SecondKmerPassTime = std::to_string(omp_get_wtime() - firstpass) + " seconds";
		printLog(SecondKmerPassTime);

		// Reliable bounds computation using estimated error rate from phred quality score
		lower = computeLower(coverage, b_pars.errorRate, b_pars.kmerSize, b_pars.minProbability);
		upper = computeUpper(coverage, b_pars.errorRate, b_pars.kmerSize, b_pars.minProbability);
		printLog(lower);
		printLog(upper);
		printLog(b_pars.errorRate);

		auto lt = countsdenovo.lock_table(); // our counting
		for (const auto &it : lt) 
		{
			if (it.second >= lower && it.second <= upper)
			{
				countsreliable_denovo.insert(it.first, kmer_id_denovo);
				++kmer_id_denovo;
			}
		}
		lt.unlock(); // unlock the table

		// Print some information about the table
		if (countsreliable_denovo.size() == 0)
		{
			std::string ErrorMessage = "BELLA terminated: 0 entries within reliable range. You may want to reduce the k-mer lenght.";
			printLog(ErrorMessage);
			exit(1);
		} 
		else 
		{
			int numReliableKmers = countsreliable_denovo.size();
			printLog(numReliableKmers);
		}

		countsdenovo.clear(); // free
		
	} // for all splits
}

#endif
