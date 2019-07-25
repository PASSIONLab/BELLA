/*********************************************************************************
 Copyright (c) 2016 Marius Erbert, Steffen Rechner

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 of the Software, and to permit persons to whom the Software is furnished to do
 so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *********************************************************************************/

#ifndef KMERHASHER_H_
#define KMERHASHER_H_

#include "ThreadBarrier.h"
#include "Bundle.h"
#include "CpuHasher.h"
#include "GpuHasher.h"
#include "KmerDistributer.h"

namespace gerbil {

/*
 * This class does the counting of the kmers. It reads SuperBundles and stores
 * counting results as KmcBundles. The class starts a number of threads that
 * use hash tables for counting.
 */
	class KmerHasher {

	private:

		// input queue for super bundles
		SyncSwapQueueSPMC<SuperBundle> *_superBundleQueue;

		// output queue for kmer count bundles
		SyncSwapQueueMPSC<KmcBundle> _kmcSyncSwapQueue;

		// threading stuff
		std::thread **_processSplitterThreads;
		std::thread *_processThread;
		uint8 _processSplitterThreadsNumber;
		uint8 _numCPUHasher;
		uint8 _numGPUHasher;
		Barrier *barrier;                                // a thread barrier

		// The kmer distributor observes throughput of the hash tables and
		// determines the ratio of kmers of a file for each hash table
		KmerDistributer *_distributor;

		// auxilliary variables
		const uint32_t _k;
		const uint_tfn _tempFilesNumber;    // number of temporary files
		TempFile *_tempFiles;                // array of temporary files
		const uint32 _thresholdMin;        // minimal occurance of kmers to be output
		uint_tfn *_tempFilesOrder;            // processing order of temporary files
		std::string _tempFolder;            // directory of temporary files ?
		const bool _norm;                    // whether to use normalized kmers

		// for synchronization
		uint32 _syncSplitterCounter;
		std::mutex _mtx;
		std::condition_variable _cv_barrier;

		// for statistic
		uint64 _kMersNumberCPU, _kMersNumberGPU;
		uint64 _uKMersNumberCPU, _uKMersNumberGPU;
		uint64 _btUKMersNumberCPU, _btUKMersNumberGPU;

		uint64 _maxKmcHashtableSize;
		uint32 _kMerBundlesNumber;

		uint64 _histogram[HISTOGRAM_SIZE];

#if false
		std::atomic<uint64> _test_kmerCounter;
		std::atomic<uint64> _testS_kmerCounter;
		std::atomic<uint64> _test_smerCounter;
#endif

		/**
		 * Spawns Splitter Threads.
		 */
		template<uint32_t K>
		void processSplit(SyncSwapQueueMPSC<cpu::KMerBundle<K>> **cpuKMerQueues,
		                  SyncSwapQueueMPSC<gpu::KMerBundle<K>> **gpuKMerQueues) {

			// spawn new threads that do the splitting
			for (uint8_t tId = 0; tId < _processSplitterThreadsNumber; ++tId) {
				_processSplitterThreads[tId] =

						new std::thread(
								[this](const uint8_t &tId, SyncSwapQueueMPSC<cpu::KMerBundle<K>> **cpuKMerQueues,
								       SyncSwapQueueMPSC<gpu::KMerBundle<K>> **gpuKMerQueues) {

									if (_norm) {
										// use normalized kmers
										processThreadSplit<K, true>(tId, cpuKMerQueues, gpuKMerQueues);
									} else {
										// do not use normalized kmers
										processThreadSplit<K, false>(tId, cpuKMerQueues, gpuKMerQueues);
									}
								}, tId, cpuKMerQueues, gpuKMerQueues);
			}
		}

		template<uint32_t K>
		void process_template() {

			_processThread =
					new std::thread([this]() {
#ifdef DEB
						printf("KmerHasher start...\n");
						StopWatch swr;
						StopWatch swt(CLOCK_THREAD_CPUTIME_ID);
						swr.start();
						swt.start();
#endif

						// create Hash Tables (CPU-side)
						cpu::HasherTask<K> cpuHasher(_numCPUHasher,
						                             _distributor,
						                             &_kmcSyncSwapQueue, _tempFiles, _thresholdMin,
						                             _maxKmcHashtableSize, _tempFolder);

						// create  Hash Tables (GPU-side)
						gpu::HasherTask<K> gpuHasher(_numGPUHasher,
						                             _distributor,
						                             &_kmcSyncSwapQueue, _tempFiles, _tempFolder, _thresholdMin);

						// Create Queues for kmer bundles
						SyncSwapQueueMPSC<cpu::KMerBundle<K>> **cpuKMerQueues;
						SyncSwapQueueMPSC<gpu::KMerBundle<K>> **gpuKMerQueues;

						cpuKMerQueues = new SyncSwapQueueMPSC<cpu::KMerBundle<K>> *[_numCPUHasher];
						gpuKMerQueues = new SyncSwapQueueMPSC<gpu::KMerBundle<K>> *[_numGPUHasher];

						const uint32_t total_threads = _numCPUHasher + _numGPUHasher;

						for (uint8 i = 0; i < _numCPUHasher; ++i) {
							// TODO: Calibrate the argument here
							cpuKMerQueues[i] = new SyncSwapQueueMPSC<cpu::KMerBundle<K>>(
									_kMerBundlesNumber / total_threads);
						}

						for (uint8 i = 0; i < _numGPUHasher; ++i) {
							// TODO: Calibrate the argument here
							gpuKMerQueues[i] = new SyncSwapQueueMPSC<gpu::KMerBundle<K>>(
									_kMerBundlesNumber / total_threads);
						}

						// start splitting
						processSplit<K>(cpuKMerQueues, gpuKMerQueues);

						// start hashing
						cpuHasher.hash(cpuKMerQueues);
						gpuHasher.hash(gpuKMerQueues);

						// wait for splitters to finalize
						joinSplitterThreads();

						// wait for queues to finalize
						for (uint8_t i = 0; i < _numCPUHasher; i++)
							cpuKMerQueues[i]->finalize();

						for (uint8_t i = 0; i < _numGPUHasher; i++)
							gpuKMerQueues[i]->finalize();

						// wait for hashtables to finalize
						cpuHasher.join();
						gpuHasher.join();

						// wait for kmc queue to finalize
						_kmcSyncSwapQueue.finalize();

						// clean up
						for (uint8 i = 0; i < _numCPUHasher; ++i)
							delete cpuKMerQueues[i];

						for (uint8 i = 0; i < _numGPUHasher; ++i)
							delete gpuKMerQueues[i];

						delete[] cpuKMerQueues;
						delete[] gpuKMerQueues;

						// statistical / debug output
						cpuHasher.printStat();
						for (uint i = 0; i < HISTOGRAM_SIZE; ++i)
							_histogram[i] = cpuHasher.getHistogramEntry(i) + gpuHasher.getHistogramEntry(i);

						_kMersNumberCPU = cpuHasher.getKMersNumber();
						_uKMersNumberCPU = cpuHasher.getUKMersNumber();
						_btUKMersNumberCPU = cpuHasher.getBtUKMersNumber();
						_kMersNumberGPU = gpuHasher.getKMersNumber();
						_uKMersNumberGPU = gpuHasher.getUKMersNumber();
						_btUKMersNumberGPU = gpuHasher.getBtUKMersNumber();
					});

		}

		/**
		 * * This macro distributes the kmer to one of the kmer bundle queues.
		 */
#define ADD_KMER_TO_BUNDLE(kmer, curTempFileId, curTempRun)                                                    \
    {                                                                                            \
        /* 																						\
		 * determine id of hash table according to current										\
		 * split ratio.																			\
		 */                                                                                        \
        const uint32_t h = _distributor->distributeKMer<K>(kmer, curTempFileId);                    \
        if(h != NULL_BUCKET_VALUE) {                                                            \
            /* assign kmer to bundle cpu or gpu */                                                    \
            if(h >= _numCPUHasher) {                                                                \
                /* add to gpu */                                                                    \
                if (!gpuKMerBundles[h-_numCPUHasher]->add(kmer)) {                                    \
                    gpuKMerBundles[h-_numCPUHasher]->setTempFileId(curTempFileId);                    \
                    gpuKMerBundles[h-_numCPUHasher]->setTempFileRun(curTempRun);                    \
                    gpuKMerQueues[h-_numCPUHasher]->swapPush(gpuKMerBundles[h-_numCPUHasher]);        \
                    gpuKMerBundles[h-_numCPUHasher]->add(kmer);                                        \
                }                                                                                    \
            }                                                                                        \
            else {                                                                                    \
                /* add to cpu */                                                                    \
                if (!cpuKMerBundles[h]->add(kmer)) {                                                \
                    cpuKMerBundles[h]->setTempFileId(curTempFileId);                                \
                    cpuKMerBundles[h]->setTempFileRun(curTempRun);                                \
                    cpuKMerQueues[h]->swapPush(cpuKMerBundles[h]);                                    \
                    cpuKMerBundles[h]->add(kmer);                                                    \
                }                                                                                    \
            }                                                                                        \
        }                                                                                             \
    }

		/**
		 * The Task of each Splitter Thread:
		 * Split SuperMers into kmers and store them into parameter queues.
		 */
		template<uint32_t K, bool NORM>
		void processThreadSplit(const uint8_t &threadId,
		                        SyncSwapQueueMPSC<cpu::KMerBundle<K>> **cpuKMerQueues,
		                        SyncSwapQueueMPSC<gpu::KMerBundle<K>> **gpuKMerQueues) {
#ifdef DEB_MESS_SUPERSPLITTER
			// do timing stuff
			StopWatch sw(CLOCK_THREAD_CPUTIME_ID);
			sw.start();
#endif

			uint64 temp_sbc(0);
			uint64 temp_kmerc(0);
			uint64 temp_kmerc2(0);

			// temporary variables
			SuperBundle *sb = new SuperBundle();// the current super bundle in use
			KMer<K> kMer, iKMer;        // the current k-mer/ inverse k-mer in use
			const KMer<K> *nKMer;                    // the current normalized k-mer
			uint16 l;                                // size of current super s-mer
			byte *b;                        // byte representation of current s-mer.
			size_t h;                // current hash value of k-mer (see distributor)
			byte nextBase;                            // next base of current s-mer
			uint_tfn curTempFileId;    // id of temp file, note: each temp file has its own id

			// Kmer bundles
			cpu::KMerBundle<K> **cpuKMerBundles =
					new cpu::KMerBundle<K> *[_numCPUHasher];
			gpu::KMerBundle<K> **gpuKMerBundles =
					new gpu::KMerBundle<K> *[_numGPUHasher];

			// init Kmer bundles
			for (uint8_t i = 0; i < _numCPUHasher; i++)
				cpuKMerBundles[i] = new cpu::KMerBundle<K>();
			for (uint8_t i = 0; i < _numGPUHasher; i++)
				gpuKMerBundles[i] = new gpu::KMerBundle<K>();

#ifdef DEB_MESS_SUPERSPLITTER
			sw.hold();	// timing
#endif

			// take next superbundle out of superbundle queue
			_superBundleQueue->swapPop(sb);

			// timing
#ifdef DEB_MESS_SUPERSPLITTER
			sw.proceed();
#endif

			uint curTempRun = 0;

#if false
			uint64 test_kmerCounter = 0;
			uint64 testS_kmerCounter = 0;
			uint64 test_smerCounter = 0;
#endif

			// for each temporary file
			for (uint rdyNmbs(0); rdyNmbs < _tempFilesNumber; ++rdyNmbs) {

				// not clear what is going on here?
				curTempFileId = _tempFilesOrder[rdyNmbs];
				if(_tempFiles[curTempFileId].isEmpty())
					continue;
				if (curTempRun + 1 < _tempFiles[curTempFileId].getNumberOfRuns()) {
					--rdyNmbs;
				}
				// get number of kmers in this file
				uint64_t kmersInFile = _tempFiles[curTempFileId].getKMersNumber();

				/* Update the distributor that controls the ratio of
				 *  kmers distributed the various hash tables. */
				if (threadId == 0)
					_distributor->updateFileInformation(curTempFileId, curTempRun,
					                                    _tempFiles[curTempFileId].getNumberOfRuns(), kmersInFile);
				this->barrier->sync();

				// Some sort of progress bar
				if (threadId == 0)
					putchar('*');

				// if super bundle is not yet empty and belongs the the current file
				if (!sb->isEmpty() && sb->tempFileId == curTempFileId && sb->tempFileRun == curTempRun) {

#ifdef DEB_MESS_SUPERSPLITTER
					// timing stuff
					sw.hold();
#endif

					// while superbundle is not yet empty and next supermer belongs
					// to the current file: read next superbundle, jump over if next binId
					while ((!sb->isEmpty() || _superBundleQueue->swapPop(sb))
					       && sb->tempFileId == curTempFileId && sb->tempFileRun == curTempRun) {

						// timing
#ifdef DEB_MESS_SUPERSPLITTER
						sw.proceed();
#endif


						// while ?
						while (sb->next(b, l)) {

#if false
							if(curTempRun == _tempFiles[curTempFileId].getNumberOfRuns() - 1) test_smerCounter++;
							if(curTempRun == _tempFiles[curTempFileId].getNumberOfRuns() - 1) testS_kmerCounter += l - K + 1;
#endif

							// extract first kmer out of supermer and add to kmer bundle
							if (NORM) {
								// normalize first
								KMer<K>::set(b, kMer, iKMer);
								nKMer = &kMer.getNormalized(iKMer);
								ADD_KMER_TO_BUNDLE((*nKMer), curTempFileId, curTempRun)

							} else {
								// just add to bundle
								kMer.set(b);
								ADD_KMER_TO_BUNDLE(kMer, curTempFileId, curTempRun);
							}

							// extract all other kmers from supermer
							for (uint i = K; i < l; ++i) {
								nextBase = (*(b + (i >> 2))
										>> (6 - ((i & 0x3) << 1))) & 0x3;
								kMer.next(nextBase);

								// add to bundle
								if (NORM) {
									iKMer.nextInv(nextBase);
									nKMer = &kMer.getNormalized(iKMer);
									ADD_KMER_TO_BUNDLE((*nKMer), curTempFileId, curTempRun);
								} else {
									ADD_KMER_TO_BUNDLE(kMer, curTempFileId, curTempRun);
								}
							}
						}

						// clean up?
						sb->clear();

						// timing
#ifdef DEB_MESS_SUPERSPLITTER
						sw.hold();
#endif
					}

					// set current file id to kmer bundles and push them to the queues
					for (uint8_t i(0); i < _numCPUHasher; ++i) {
						if (!cpuKMerBundles[i]->isEmpty()) {
							cpuKMerBundles[i]->setTempFileId(curTempFileId);
							cpuKMerBundles[i]->setTempFileRun(curTempRun);
							cpuKMerQueues[i]->swapPush(cpuKMerBundles[i]);
						}
					}

					// for gpu: set current file id to kmer bundles and push them to the queues
					for (uint8_t i(0); i < _numGPUHasher; ++i) {
						if (!gpuKMerBundles[i]->isEmpty()) {
							gpuKMerBundles[i]->setTempFileId(curTempFileId);
							gpuKMerBundles[i]->setTempFileRun(curTempRun);
							gpuKMerQueues[i]->swapPush(gpuKMerBundles[i]);
						}
					}

					if (!sb->isEmpty()) {
						curTempRun = sb->tempFileRun;
						//std::cout << "next run: " << curTempRun << std::endl;
						//std::cout << "next tfid: " << sb->tempFileId << std::endl;
					}

					// timing
#ifdef DEB_MESS_SUPERSPLITTER
					sw.proceed();
#endif
				}

				//MEMORY_BARRIER
				memoryBarrier();

				this->barrier->sync();
			}

#if false
			_test_kmerCounter += test_kmerCounter;
			_testS_kmerCounter += testS_kmerCounter;
			_test_smerCounter += test_smerCounter;
#endif

			// clean up
			delete sb;

			// free kmer bundles
			for (uint8_t i = 0; i < _numCPUHasher; i++)
				delete cpuKMerBundles[i];
			for (uint8_t i = 0; i < _numGPUHasher; i++)
				delete gpuKMerBundles[i];

			delete[] cpuKMerBundles;
			delete[] gpuKMerBundles;

// timing
#ifdef DEB_MESS_SUPERSPLITTER
			sw.stop();
			printf("superplitter[%2u]: %.3f s\n", threadId, sw.get_s());
#endif
		}

		/**
		 * Wait for all splitter threads to finalize.
		 */
		void joinSplitterThreads() {
			for (uint8_t i = 0; i < _processSplitterThreadsNumber; i++) {
				_processSplitterThreads[i]->join();
				delete _processSplitterThreads[i];
			}
		}

	public:

		/**
		 * Constructor.
		 */
		KmerHasher(const uint32_t &k, const uint32 &kmcBundlesNumber,
		           SyncSwapQueueSPMC<SuperBundle> *const superBundleQueue,
		           const uint8 &processSplitterThreadsNumber,
		           const uint8 &processHasherThreadsNumber,
		           const uint8 &processGPUHasherThreadsNumber, TempFile *tempFiles,
		           const uint_tfn &tempFilesNumber, const uint32 &thresholdMin,
		           const bool &norm, std::string pTempFolder,
		           const uint64 &maxKmcHashtableSize, const uint32 &kMerBundlesNumber,
		           uint_tfn *tempFilesOrder,
		           KmerDistributer *distributor
		) :
				_kmcSyncSwapQueue(kmcBundlesNumber),
				_processSplitterThreadsNumber(processSplitterThreadsNumber),
				_numCPUHasher(processHasherThreadsNumber),
				_numGPUHasher(processGPUHasherThreadsNumber),
				barrier(nullptr), _k(k), _tempFilesNumber(tempFilesNumber),
				_thresholdMin(thresholdMin), _norm(norm), _tempFolder(pTempFolder),
				_superBundleQueue(superBundleQueue), _tempFiles(
				tempFiles), _processThread(NULL), _kMersNumberCPU(0), _kMersNumberGPU(
				0), _uKMersNumberCPU(0), _uKMersNumberGPU(0), _btUKMersNumberCPU(
				0), _btUKMersNumberGPU(0), _maxKmcHashtableSize(
				maxKmcHashtableSize), _syncSplitterCounter(0), _kMerBundlesNumber(
				kMerBundlesNumber), _tempFilesOrder(tempFilesOrder), _distributor(distributor) {

			// create array of threads
			_processSplitterThreads =
					new std::thread *[_processSplitterThreadsNumber];

			barrier = new Barrier(_processSplitterThreadsNumber);

#if false
			_test_kmerCounter = 0;
			_testS_kmerCounter = 0;
			_test_smerCounter = 0;
#endif
		}

		/**
		 * Clean up.
		 */
		~KmerHasher() {
			delete[] _processSplitterThreads;
			delete barrier;

#if false
			printf("kmerHasher kmers: %lu\n", _test_kmerCounter.load());
			printf("kmerHasher S kmers: %lu\n", _testS_kmerCounter.load());
			printf("kmerHasher smers: %lu\n", _test_smerCounter.load());
#endif
		}

		void saveHistogram() {
			FILE *file;
			file = fopen((_tempFolder + "histogram.csv").c_str(), "wb");
			fprintf(file, "counter; number of uk-mers\n");
			for (uint i = 1; i < HISTOGRAM_SIZE; ++i)
				fprintf(file, "  %3u; %9lu\n", i, _histogram[i]);
			fprintf(file, ">=%3u; %9lu\n", HISTOGRAM_SIZE, _histogram[0]);
			fclose(file);
		}

		/** Main working procedure for this class. */
		void process() {

#define C_PROC(x) case x: process_template<x>(); break

			// decide which template specialization to use
			switch (_k) {

					//LOOP512(MAX_KMER_SIZE, C_PROC);
					LOOP128(MAX_KMER_SIZE, C_PROC);
					//C_PROC(28);
					//C_PROC(40);
					//C_PROC(56);
					//C_PROC(65);
				default:
					throw std::runtime_error(
							std::string("Gerbil Error: Unsupported k"));
			}
		}

		void join() {
			_processThread->join();
			delete _processThread;
			IF_DEB(printf("all KmerHashers are rdy...\n"));
		}

		void print() {
			printf("kmers (CPU)     : %12lu\n", _kMersNumberCPU);
			printf("kmers (GPU)     : %12lu\n", _kMersNumberGPU);
			printf("ukmers (CPU)    : %12lu\n", _uKMersNumberCPU);
			printf("ukmers (GPU)    : %12lu\n", _uKMersNumberGPU);
			printf("below th (CPU)  : %12lu\n", _btUKMersNumberCPU);
			printf("below th (GPU)  : %12lu\n", _btUKMersNumberGPU);
		}

		SyncSwapQueueMPSC<KmcBundle> *getKmcSyncSwapQueue() {
			return &_kmcSyncSwapQueue;
		}
	};

}

#endif /* KMERHASHER_H_ */
