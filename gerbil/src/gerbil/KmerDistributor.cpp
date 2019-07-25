/*
 * KmerDistributer.cpp
 *
 *  Created on: Jan 22, 2016
 *      Author: rechner
 */

#include "../../include/gerbil/KmerDistributer.h"
#include "../../include/gerbil/config.h"
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>
#include <map>
#include <algorithm>
#include <sstream>

namespace gerbil {

	KmerDistributer::KmerDistributer(const uint32_t numCPUHasherThreads,
	                                 const uint32_t numGPUHasherThreads, const uint32_t maxNumFiles) :
			numThreadsCPU(numCPUHasherThreads), numThreadsGPU(numGPUHasherThreads),
			_kmerNumber(0), _ukmerNumber(0), _ukmerRatio(0.0), _totalCapacity(0) {

		const uint32_t numThreads = numCPUHasherThreads + numGPUHasherThreads;

		prefixSum = new double[numThreads + 1];
		tmp = new double[numThreads];
		capacity = new uint64_t[numThreads];

		throughput = new double[numThreads];
		memset(throughput, 0, numThreads * sizeof(double)); // fill throughput array with initial zeros

		// init all capacities with infinity
		for (int i = 0; i < numThreads; i++)
			capacity[i] = ULONG_MAX;

		// initialize file specific variables
		ratio.resize(maxNumFiles);
		buckets.resize(maxNumFiles);
		lock.resize(maxNumFiles);
		for (int i = 0; i < maxNumFiles; i++) {
			ratio[i] = new double[numThreads];
			buckets[i] = new uint32_t[BUCKETSIZE];
			lock[i] = 1;                // set lock status of each file to locked
			memset(ratio[i], 0, numThreads * sizeof(double));
		}
	}

	KmerDistributer::~KmerDistributer() {
		delete[] prefixSum;
		delete[] tmp;
		delete[] capacity;
		delete[] throughput;
		for (int i = 0; i < buckets.size(); i++) {
			delete[] buckets[i];
			delete[] ratio[i];
		}
	}

	void KmerDistributer::updateThroughput(const bool gpu, const uint32_t tId,
	                                       const double t) {

		if (gpu) {
			//printf("gpu hasher %u: update throughput to %f kmers/ms\n", tId, t);
			throughput[numThreadsCPU + tId] = t;
		} else {
			//printf("cpu hasher %u: update throughput to %f kmers/ms\n", tId, t);
			throughput[tId] = t;
		}
	}

	void KmerDistributer::updateCapacity(const bool gpu, const uint32_t tId,
	                                     const uint64_t cap) {
		if (gpu) {
			//printf("gpu hasher %u: update capacity to %u kmer\n", tId, cap);
			capacity[numThreadsCPU + tId] = cap;
		} else {
			//printf("cpu hasher %u: update capacity to %u kmer\n", tId, cap);
			//capacity[tId] = cap;
		}
		_totalCapacity += cap;
		//std::cout << "cap: " << cap  << "   " << _totalCapacity << std::endl;
		//printf("updateCapacity: cap= %lu   _totalCapacity= %lu\n", cap, _totalCapacity.load());

	}


	/*void KmerDistributer::updateUKmerRatio(const bool gpu, const uint32_t tId, const uint binkmernumber, const uint binukmernumber) {
		if(!tId && !gpu) {
			_ukmerNumber += binukmernumber;
			_kmerNumber += binkmernumber;
			_ukmerRatio = _kmerNumber ? ((double)_ukmerNumber) / _kmerNumber : 1;
		}
	}*/

	double KmerDistributer::getSplitRatio(const bool gpu, const uint32_t tId,
	                                     const uint_tfn curTempFileId) const {

		//return 1;

		// wait until updates are complete
		//while (lock[curTempFileId])
		//	usleep(10);

		return gpu ?
		       ratio[curTempFileId][numThreadsCPU + tId] :
		       ratio[curTempFileId][tId];
	}

	void KmerDistributer::updateFileInformation(const uint32_t curFileId, const uint32_t curRun, const uint32_t maxRunsPerFile,
	                                            const uint64_t kmersInFile) {
		//printf(">>>>>>>>>>>>>updateFileInformation!\n");

		// try to lock data structure
		//if (CAS(&lock[curFileId], 1, 2) != 1)
		//	return;

		const double eps = 1e-6;
		const uint32_t numThreads = numThreadsCPU + numThreadsGPU;

		// compute sum of all throughput
		double sum = 0;
		for (int i = 0; i < numThreads; i++) {
			//printf("thread %u has throughput %f\n", i, throughput[i]);
			sum += throughput[i];

			// break if any thread has not yet provieded throughput information
			if (fabs(throughput[i]) < eps) {
				sum = 0;
				break;
			}
		}

		/**
		 * Update the ratio of kmers from latest throughput and capacity
		 */

		// if some thread has not provided throughput information
		if (fabs(sum) < eps) {
			// distribute kmers equally over all threads
			for (int i = 0; i < numThreads; i++) {
				ratio[curFileId][i] = 1.0f / numThreads;
				//printf("thread %u becomes ratio of %f\n", i, ratio[curFileId][i]);
			}
		}
			// if all threads have provided throughput information
		else {

			// define the new ratio of kmers for each hasher thread
			for (int i = 0; i < numThreads; i++) {
				/* Version 1: Each Ratio is defined as average between throughput ratio and uniform distribution */
				//cpuRatio[i] = 0.5 * cpu_throughput[i] / sum + 0.5 / (numCPUHasherThreads + numGPUHasherThreads);
				/* Version 2: Each Ratio is defined as average between troughput ratio and old ratio */
				ratio[curFileId][i] = 0.5 * ratio[lastFileId][i]
				                      + 0.5 * throughput[i] / sum;

				/* Version 3: Each Ratio is defined as pure throughput ratio */
				//ratio[curFileId][i] = throughput[lastFileId][i] / sum;
				//printf("thread %u becomes ratio of %f\n", i, ratio[curFileId][i]);
			}
		}

		// Check whether any capacity constraint will be violated.
		int numPositiveCapacity = 0;
		for (int i = 0; i < numThreads; i++) {
			double maxRatio = (double) capacity[i] / (double) kmersInFile;
			tmp[i] = maxRatio - ratio[curFileId][i];
			if (tmp[i] > 0)
				numPositiveCapacity++;
			//printf("thread %u has diff of %f\n", i, tmp[i]);
		}

		for (int i = 0; i < numThreads; i++) {

			// while ratio constraint of thread i is violated and there are threads that could compensate
			int rnd = 0;
			while (rnd < 100 && tmp[i] < -eps && numPositiveCapacity > 0) {

				rnd++;

				//printf("while (tmp[i] < 0 && numPositiveCapacity > 0)\n");

				// try to redistribute an amount of kmers to other hasher threads
				/*printf(
						"capacity of thread %u is violated! redistribute a relative amount of %f kmers!\n",
						i, -tmp[i]);*/

				// try to redistribute uniformly to all threads with positive capacity
				double y = -tmp[i] / (double) numPositiveCapacity;

				for (int k = 0; tmp[i] < -eps && k < numThreads; k++) {

					// if thread k has still open capacity
					if (tmp[k] > eps) {

						double x = std::min(y, tmp[k]);
						ratio[curFileId][i] -= x;// ratio of thread i becomes smaller by x
						ratio[curFileId][k] += x;// ratio of thread k becomes larger by x
						tmp[i] += x;
						tmp[k] -= x;
						//printf("give %f to thread %u\n", x, k);
						if (tmp[k] <= eps)
							numPositiveCapacity--;
					}
				}
			}

			// if difference is still smaller than zero: exception!
			if (numPositiveCapacity <= 0 && tmp[i] < -eps)
				throw std::runtime_error(
						"gerbil::exception: system has not enough memory for this number of files. Try to increase the number of temporary files.\n");
		}

		/**
		 *  Update prefix sums
		 */
		int j = 0;
		double x = 0;
		for (int i = 0; i < numThreads; i++) {
			prefixSum[j] = x;
			x += ratio[curFileId][i];
			j++;
			//	printf("thread %u becomes final ratio of %f\n", i, ratio[curFileId][i]);
		}
		prefixSum[numThreads] = 1.0;

		/**
		 * Update Bucket table for the current file.
		 */
		uint32_t* bucket = buckets[curFileId];
		for(size_t i = 0; i < BUCKETSIZE; ++i)
			bucket[i] = NULL_BUCKET_VALUE;

		uint tempBucketSize = BUCKETSIZE / maxRunsPerFile;
		uint offset = curRun * tempBucketSize;

		if(curRun + 1 == maxRunsPerFile)
			tempBucketSize = BUCKETSIZE - offset;

		uint32_t threadId = 0;
		int i = 0;
		// for each hasher thread
		while (threadId < numThreads) {

			// for each bucket entry with corresponds to current hasher thread
			//std::cout << "tempBucketSize: " << tempBucketSize << "   " << (int)(prefixSum[bucketId + 1] * tempBucketSize) << std::endl;
			while (offset + i < BUCKETSIZE && i < prefixSum[threadId + 1] * tempBucketSize) {
				bucket[offset + i] = threadId;
				i++;
			}
			threadId++;
		}

#if false
		std::stringstream s;
		s << (int) curFileId << "  " << (int) curRun << "  [";
		for(size_t i = 0; i < BUCKETSIZE; ++i) {
			if(bucket[i] == NULL_BUCKET_VALUE)
				s << "N";
			else
				s << (int) bucket[i];
			s << ", ";
		}
		s << "]\n";
		std::cout << s.str() << std::flush;
#endif
		// unlock current file
		//lock[curFileId] = 0;

		lastFileId = curFileId;
	}

} /* namespace gerbil */
