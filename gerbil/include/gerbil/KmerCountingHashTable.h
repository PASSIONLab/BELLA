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

#ifndef KMERCOUNTINGHASHTABLE_H_
#define KMERCOUNTINGHASHTABLE_H_

#include "../cuda_ds/CountingHashTable.h"
#include "KMer.h"
#include "SyncQueue.h"
#include "Bundle.h"

namespace gerbil {
	namespace gpu {

		namespace cd = ::cuda_ds;
		namespace cdi = ::cuda_ds::internal;

/**
 * Extend the cuda_ds::CountingHashTable to specialize the counting of KMers.
 * Just wrapped the constructor, the add method and added an additional
 * extracting method.
 *
 * @param K K
 * @param intsPerKMer The number of ints per KMer
 */
		template<uint32_t K>
		class KmerCountingHashTable : public cd::CountingHashTable<KMer < K>,
		GPU_KMER_BUNDLE_DATA_SIZE> {

		private:
		typedef cdi::KeyValuePair<cdi::intsPerKey<KMer < K>>()>
		_KeyValuePair;

		SyncSwapQueueMPSC <KmcBundle> *kmcQueue;    // output queue
		uint32_t threshold;                        // minimal number of kmer occurence
		uint64_t kMersNumber;                    // number of inserted kmers
		uint64_t uKMersNumber;                    // number of inserted unique kmers
		uint64_t btUKMersNumber;                // number of inserted unique kmers below threshold

		char *rawbuffer;                        // buffer for copy jobs

		// buffers needed for emergency evacuation of the noSuccessArea of the Hash table
		FailureBuffer <K> *failureBuffer;
		FailureBuffer <K> *tempBuffer;
		KMerBundle <K> *_kmb;               // working kmer bundle
		KmcBundle *curKmcBundle;            // working kmc bundle

		uint64 _histogram[HISTOGRAM_SIZE];

		public:

		/**
		 * Main Constructor.
		 * @param devID the cuda device id on which the table should be allocated.
		 * @param kmcQueue The output queue where kmer counts are stored.
		 * @param threshold The minimal number of counts to be output.
		 */
		KmerCountingHashTable(const uint32_t devID,
		                      SyncSwapQueueMPSC <KmcBundle> *kmcQueue, uint32_t threshold,
		                      std::string tempPath) :
				cd::CountingHashTable<KMer < K>,GPU_KMER_BUNDLE_DATA_SIZE>(devID),

		kmcQueue (
		kmcQueue), threshold(threshold), kMersNumber(0), uKMersNumber(
				0), btUKMersNumber(0), rawbuffer(nullptr), failureBuffer(
				nullptr), tempBuffer(nullptr), _kmb(nullptr), curKmcBundle(
				nullptr) {

			// allocate host memory buffer
			rawbuffer = (char *) malloc(GPU_COPY_BUFFER_SIZE);
			failureBuffer = new FailureBuffer<K>(
					FAILUREBUFFER_KMER_BUNDLES_NUMBER_PER_THREAD, tempPath, 1024 + devID,
					0);
			tempBuffer = new FailureBuffer<K>(
					FAILUREBUFFER_KMER_BUNDLES_NUMBER_PER_THREAD, tempPath, 1024 + devID,
					1);
			_kmb = new KMerBundle<K>();
			curKmcBundle = new KmcBundle();

			// initialize histogram
			uint64 histogram[HISTOGRAM_SIZE];
			for (size_t i(0); i < HISTOGRAM_SIZE; ++i)
				histogram[i] = 0;
		}

		~KmerCountingHashTable() {
			free(rawbuffer);
			delete failureBuffer;
			delete tempBuffer;
			delete _kmb;
			delete curKmcBundle;
		}

		uint64 getHistogramEntry(const uint i) {
			return _histogram[i];
		}

		/**
		 * Get number of kmers inserted until last call of extract().
		 */
		uint64_t getKMersNumber() const {
			return kMersNumber;
		}

		/**
		 * Get number of distinct kmers inserted until last call of extract().
		 */
		uint64_t getUKMersNumber() const {
			return uKMersNumber;
		}

		/**
		 * Get number of distinct kmers below the threshold until last call of extract().
		 */
		uint64_t getUKMersNumberBelowThreshold() const {
			return btUKMersNumber;
		}

		/**
		 * Add a bundle of kmers to the table.
		 * @param kmb Pointer to the kmer bundle to be inserted.
		 * @return True, if the kmers in this bundle could be stored successfully.
		 * False, if an emergency extract is neccessary.
		 */
		void addBundle(gpu::KMerBundle<K> *kmb) {

			const KMer<K> *data = (const KMer<K> *) kmb->getData();
			uint64_t kmb_count = kmb->count();

			// try to insert the bundle into the hash table
			bool success = this->insert(data, kmb_count);

			if (!success) {

				printf("EMERGENCY!\n");
				emergencyExtract();

				// try to add again
				this->insert((const KMer<K> *) kmb->getData(), kmb->count());
			}

			kmb->clear();
		}

		/**
		 * Extract the kmers with their counts and push them into the kmc bundle queue.
		 */
		void extract() {

			uint64_t extractedUKmers = uKMersNumber;

			// define pointers and sizes
			KMer<K> *kmers = (KMer<K> *) rawbuffer;                         // kmer-based view on the buffer
			_KeyValuePair *kmerCounts = (_KeyValuePair * )rawbuffer;        // entry-based view on the buffer
			const uint32_t kmersPerBuffer = GPU_COPY_BUFFER_SIZE / sizeof(KMer<K>);     // number of kmers that fit in the buffer
			const uint64_t sizeOfEntry = sizeof(_KeyValuePair);
			const uint32_t entriesPerBuffer = GPU_COPY_BUFFER_SIZE / sizeOfEntry;

			// temporary working variables
			KMer<K> *curKey;
			_KeyValuePair * curKeyValuePair;
			uint32_t curCount;

			bool done = false;
			while(!done) {    // run as long as some failures exist

				// set device id
				cudaSetDevice(this->devID);

				// wait for current threads to finish their work
				this->join();

				/**
				 * 1. Extract kmers and counts in hash table.
				 * Copy the whole content of the hash table back to main memory.
				 */

				// compress table entries
				//const uint64_t numEntries = this->compressEntries();
				//printf("before compression: num entries=%lu\n", this->getNumEntries());
				//printf("after compression: numEntries=%lu\n", numEntries);
				const uint64_t numEntries = this->getNumEntries();

				// Load table to host memory in small chunks of size GPU_COPY_BUFFER_SIZE
				for (uint32_t i = 0; i * entriesPerBuffer < numEntries;  i++) {

					const uint64_t itemsToFetch =
							(i + 1) * entriesPerBuffer <= numEntries ?
							entriesPerBuffer :
							numEntries % entriesPerBuffer;

					// copy chunk of data
					cudaMemcpy(rawbuffer, &this->table[i * entriesPerBuffer],
					           itemsToFetch * sizeOfEntry, cudaMemcpyDeviceToHost);

					// copy kmer counts to result vector
					for (uint32_t j = 0; j < itemsToFetch; j++) {

						curKeyValuePair = (_KeyValuePair * ) & kmerCounts[j];
						curKey = (KMer<K> *) &curKeyValuePair->getKey();
						curCount = curKeyValuePair->getCount();

						if (!curKey->isEmpty() && curCount > 0) {
							++_histogram[curCount < HISTOGRAM_SIZE ? curCount : 0];
							if (curCount >= threshold) {

								// copy kmer count to output queue
								if (!curKmcBundle->add<K>(*curKey, curCount)) {
									kmcQueue->swapPush(curKmcBundle);
									curKmcBundle->add<K>(*curKey, curCount);
								}
								curKey->clear();
							} else {
								++btUKMersNumber;
							}

							++uKMersNumber;
							kMersNumber += curCount;
						}
					}
				}

				/**
				 * 2. Extract kmers in noSuccess area.
				 * First, try to sort them in gpu memory and copy them back to main memory.
				 * If something goes wrong with sorting, copy the whole thing back,
				 * reset the table and insert them again.
				 */

				// Determine if there are any unsuccessfully inserted kmers
				uint64_t numNoSuccess;
				cudaMemcpy(&numNoSuccess, this->numNoSuccessPtr, sizeof(uint64_t), cudaMemcpyDeviceToHost);

				//printf("numNoSuccess=%lu\n", numNoSuccess);

				// if there are some kmer in the noSuccess area
				if (numNoSuccess > 0) {

					// the flag allgood is true, if entries in noSuccess area are free to copy
					// them back to main memory
					bool allgood = failureBuffer->isEmpty();

					// problem occur when sorting to many keys (kernel is being killed)
					if(numNoSuccess > 50*1000*1000) {
						allgood = false;
					}
					else {
						// try to sort keys and copy them back to kmc bundle
						try {
							cdi::sortKeys<cdi::intsPerKey<KMer<K>>()>(this->noSuccessArea, numNoSuccess);
						} catch (...) {
							// clear cuda Error stack
							if(cudaGetLastError() != cudaSuccess) {
								allgood = false;  // if kmers could not be sorted by any reason
							}
							// if kmers could not be sorted by any reason
							allgood = false;
						}
					}

					if(allgood) {   // area is clean and could be sorted

						// working variables
						uint32_t curCount;

						// load first buffer
						uint64_t itemsToFetch =
								kmersPerBuffer <= numNoSuccess ?
								kmersPerBuffer : numNoSuccess % kmersPerBuffer;

						cudaMemcpy(rawbuffer, this->noSuccessArea,
						           itemsToFetch * sizeof(KMer<K>), cudaMemcpyDeviceToHost);

						// compare with last seen key
						KMer<K> lastKey = kmers[0];
						curCount = 1;

						// compress and add to result vector
						for (uint64_t j = 1; j < itemsToFetch; j++) {

							if (kmers[j] != lastKey) {

								if (!lastKey.isEmpty()) {
									++_histogram[curCount < HISTOGRAM_SIZE ? curCount : 0];
									if (curCount >= threshold) {

										// copy kmer count to output queue
										if (!curKmcBundle->add<K>(lastKey, curCount)) {
											kmcQueue->swapPush(curKmcBundle);
											curKmcBundle->add<K>(lastKey, curCount);
										}
									} else {
										++btUKMersNumber;
									}

									++uKMersNumber;
									kMersNumber += curCount;
								}

								lastKey = kmers[j];
								curCount = 1;
							} else
								curCount++;
						}

						// load other buffers (if any)
						for (uint64_t i = 1; i * kmersPerBuffer < numNoSuccess; i++) {

							itemsToFetch =
									(i + 1) * kmersPerBuffer <= numNoSuccess ?
									kmersPerBuffer : numNoSuccess % kmersPerBuffer;

							cudaMemcpy(rawbuffer, &this->noSuccessArea[i * kmersPerBuffer],
							           itemsToFetch * sizeof(KMer<K>),
							           cudaMemcpyDeviceToHost);

							// compress and add to result vector
							for (uint64_t j = 0; j < itemsToFetch; j++) {
								if (kmers[j] != lastKey) {

									// put out kmer count
									if (!lastKey.isEmpty()) {
										++_histogram[curCount < HISTOGRAM_SIZE ? curCount : 0];
										if (curCount >= threshold) {

											// copy kmer count to output queue
											if (!curKmcBundle->add<K>(lastKey, curCount)) {
												kmcQueue->swapPush(curKmcBundle);
												curKmcBundle->add<K>(lastKey, curCount);
											}
										} else {
											btUKMersNumber++;
										}

										uKMersNumber++;
										kMersNumber += curCount;
									}

									lastKey = kmers[j];
									curCount = 1;
								} else
									curCount++;
							}
						}

						// insert last seen kmer
						if (!lastKey.isEmpty()) {
							++_histogram[curCount < HISTOGRAM_SIZE ? curCount : 0];
							if (curCount >= threshold) {
								// copy kmer count to output queue
								if (!curKmcBundle->add<K>(lastKey, curCount)) {
									kmcQueue->swapPush(curKmcBundle);
									curKmcBundle->add<K>(lastKey, curCount);
								}
							} else {
								btUKMersNumber++;
							}

							uKMersNumber++;
							kMersNumber += curCount;
						}

					}
					else {

						// if kmers could not be sorted or a failure occurred during insertion

						// copy the whole noSuccess area back to main memory
						// and store the kmers in failure buffers

						// for all chunks
						for (uint64_t i = 0; i * kmersPerBuffer < numNoSuccess; i++) {

							const uint64_t itemsToFetch =
									(i + 1) * kmersPerBuffer <= numNoSuccess ?
									kmersPerBuffer : numNoSuccess % kmersPerBuffer;

							// copy back to host memory
							cudaMemcpy(kmers, &this->noSuccessArea[i * kmersPerBuffer],
							           itemsToFetch * sizeof(KMer<K>),
							           cudaMemcpyDeviceToHost);

							// store kmers in failure buffer
							for (uint32_t j = 0; j < itemsToFetch; j++) {
								failureBuffer->addKMer(kmers[j]);
							}
						}
					}
				}

				// submit last kmc bundle
				if (!curKmcBundle->isEmpty()) {
					kmcQueue->swapPush(curKmcBundle);
				}

				/**
				 * 3. Handle failures.
				 */

				if (!failureBuffer->isEmpty()) {

					//printf("extractAndClear: handle failures\n");

					// clear table
					const uint64_t size = (uint64_t) (failureBuffer->getAmount() / FILL_GPU);
					this->init(size);

					// the original failureBuffer may be filled again
					std::swap(tempBuffer, failureBuffer);

					// for all kmer bundles in failure buffer
					while (tempBuffer->getNextKMerBundle(_kmb)) {
						// add bundle to the table
						addBundle(_kmb);
					}

					// all kmer are read from failure buffer
					tempBuffer->clear();
				}
				else {
					// table has been completely extracted
					done = true;
				}
			}
		}

		void emergencyExtract() {

			// set device id
			cudaSetDevice(this->devID);

			// wait for all threads to finish their work
			this->join();

			// entry-based view on the buffer
			KMer<K> *kmers = (KMer<K> *) rawbuffer;

			// Determine if there are any unsuccessfully inserted kmers
			uint64_t numNoSuccess;
			cudaMemcpy(&numNoSuccess, this->numNoSuccessPtr, sizeof(uint64_t),
			           cudaMemcpyDeviceToHost);

			//printf("copy %u kmers back to host\n", numNoSuccess);

			// if there are some
			if (numNoSuccess > 0) {

				// Load kmers (in small chunks) to host memory
				const uint64_t kmersPerBuffer = GPU_COPY_BUFFER_SIZE / sizeof(KMer<K>);

				// for all chunks
				for (uint32_t i = 0; i * kmersPerBuffer < numNoSuccess; i++) {

					const uint64_t itemsToFetch =
							(i + 1) * kmersPerBuffer <= numNoSuccess ?
							kmersPerBuffer : numNoSuccess % kmersPerBuffer;

					// copy back to host memory
					cudaMemcpy(kmers, &this->noSuccessArea[i * kmersPerBuffer],
					           itemsToFetch * sizeof(KMer<K>),
					           cudaMemcpyDeviceToHost);

					// store kmers in failure buffer
					for (uint32_t j = 0; j < itemsToFetch; j++) {
						failureBuffer->addKMer(kmers[j]);
					}
				}

				// reset the no Success counter
				cudaError_t err = cudaMemset(this->numNoSuccessPtr, 0, sizeof(uint64_t));
				if (err != cudaSuccess) {
					throw std::runtime_error(
							std::string("cuda_ds::emergencyExtract::exception: ")
							+ std::string(cudaGetErrorName(err)) + " "
							+ std::string(cudaGetErrorString(err)));
				}
			}
		}

#if false
		/**
		 * Extract, compress and output the entries stored at the noSuccessArea
		 * and push them as kmer counts into the kmc bundle queue.
		 */
		void extractNoSuccessArea() {

			// set device id
			cudaSetDevice(this->devID);

			// entry-based view on the buffer
			KMer<K>* kmers = (KMer<K>*) rawbuffer;

			// working variables
			uint32_t curCount;

			// Determine if there are any unsuccessfully inserted kmers
			uint64_t numNoSuccess;
			cudaMemcpy(&numNoSuccess, this->numNoSuccessPtr, sizeof(uint64_t), cudaMemcpyDeviceToHost);

			printf("numNoSuccess=%lu\n", numNoSuccess);

			// if there are some
			if (numNoSuccess > 0) {

				try {
					// sort area of shame
					cdi::sortKeys<cdi::intsPerKey<KMer<K>>()>(this->noSuccessArea, numNoSuccess);
				} catch(std::bad_alloc) {
					emergencyExtract();
					return;
				}

				// Load kmers (in small chunks)
				const uint32_t kmersPerBuffer = GPU_COPY_BUFFER_SIZE
						/ sizeof(KMer<K> );

				// load first buffer
				uint64_t itemsToFetch =
						kmersPerBuffer <= numNoSuccess ?
								kmersPerBuffer : numNoSuccess % kmersPerBuffer;

				cudaMemcpy(rawbuffer, this->noSuccessArea,
						itemsToFetch * sizeof(KMer<K> ), cudaMemcpyDeviceToHost);

				// compare with last seen key
				KMer<K> lastKey = kmers[0];
				curCount = 1;

				// compress and add to result vector
				for (uint64_t j = 1; j < itemsToFetch; j++) {

					if (kmers[j] != lastKey) {

						if (!lastKey.isEmpty()) {
							++_histogram[curCount < HISTOGRAM_SIZE ? curCount : 0];
							if (curCount >= threshold) {

								// copy kmer count to output queue
								if (!curKmcBundle->add<K>(lastKey, curCount)) {
									kmcQueue->swapPush(curKmcBundle);
									curKmcBundle->add<K>(lastKey, curCount);
								}
							} else {
								++btUKMersNumber;
							}

							++uKMersNumber;
							kMersNumber += curCount;
						}

						lastKey = kmers[j];
						curCount = 1;
					} else
						curCount++;
				}

				// load other buffers (if any)
				for (uint64_t i = 1; i * kmersPerBuffer < numNoSuccess; i++) {

					itemsToFetch =
							(i + 1) * kmersPerBuffer <= numNoSuccess ?
									kmersPerBuffer : numNoSuccess % kmersPerBuffer;

					cudaMemcpy(rawbuffer, &this->noSuccessArea[i * kmersPerBuffer],
							itemsToFetch * sizeof(KMer<K> ),
							cudaMemcpyDeviceToHost);

					// compress and add to result vector
					for (uint64_t j = 0; j < itemsToFetch; j++) {
						if (kmers[j] != lastKey) {

							// put out kmer count
							if (!lastKey.isEmpty()) {
								++_histogram[curCount < HISTOGRAM_SIZE ? curCount : 0];
								if (curCount >= threshold) {

									// copy kmer count to output queue
									if (!curKmcBundle->add<K>(lastKey, curCount)) {
										kmcQueue->swapPush(curKmcBundle);
										curKmcBundle->add<K>(lastKey, curCount);
									}
								} else {
									btUKMersNumber++;
								}

								uKMersNumber++;
								kMersNumber += curCount;
							}

							lastKey = kmers[j];
							curCount = 1;
						} else
							curCount++;
					}
				}

				// insert last seen kmer
				if (!lastKey.isEmpty()) {
					++_histogram[curCount < HISTOGRAM_SIZE ? curCount : 0];
					if (curCount >= threshold) {
						// copy kmer count to output queue
						if (!curKmcBundle->add<K>(lastKey, curCount)) {
							kmcQueue->swapPush(curKmcBundle);
							curKmcBundle->add<K>(lastKey, curCount);
						}
					} else {
						btUKMersNumber++;
					}

					uKMersNumber++;
					kMersNumber += curCount;
				}
			}
		}

		/**
		 * Extract kmer counts from the hash table and insert them
		 * as kmc bundles into the output queue. In the end, the
		 * table is empty again and has same size.
		 */
		void extractAndClear() {

			printf("extractAndClear\n");

			// wait for threads to success
			cudaError_t err = cudaDeviceSynchronize();
			if (err != cudaSuccess) {
				throw std::runtime_error(
						std::string("gerbil::exception: ")
						+ std::string(cudaGetErrorName(err)) + " "
						+ std::string(cudaGetErrorString(err)));
			}

			extractTableEntries();
			extractNoSuccessArea();

			printf("extractAndClear: clear the table\n");

			// clear the table
			this->clear();

			// handle failures as long as some exist
			while (!failureBuffer->isEmpty()) {

				printf("extractAndClear: handle failures\n");

				// the original failureBuffer may be filled again
				std::swap(tempBuffer, failureBuffer);

				// for all kmer bundles in failure buffer
				while (tempBuffer->getNextKMerBundle(_kmb)) {
					// add bundle to the table
					addBundle(_kmb);
				}

				printf("extractAndClear: all bundles added\n");

				// wait for threads to success
				cudaError_t err = cudaDeviceSynchronize();
				if (err != cudaSuccess) {
					throw std::runtime_error(
							std::string("gerbil::exception: ")
							+ std::string(cudaGetErrorName(err)) + " "
							+ std::string(cudaGetErrorString(err)));
				}

				// extract kmer counts
				extractTableEntries();
				extractNoSuccessArea();

				this->clear();

				// all kmer are read from failure buffer
				tempBuffer->clear();
			}

			// submit last kmc bundle
			if (!curKmcBundle->isEmpty()) {
				kmcQueue->swapPush(curKmcBundle);
			}

			printf("extractAndClear done\n");
		}
#endif
	};

}
}

#endif /* KMERCOUNTINGHASHTABLE_H_ */
