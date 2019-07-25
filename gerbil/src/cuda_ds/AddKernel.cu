#include "../../include/cuda_ds/CountingHashTable/_CountingHashTable.h"

#include <iostream>
#include <stdexcept>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>

namespace cuda_ds {
namespace internal {

//#define DEBUG

// constants
const uint32_t MAXTRIALS = 30;
const uint32_t STATUS_LOCKED = 1;
const uint32_t STATUS_FREE = 2;
const uint32_t STATUS_OCCUPIED = 3;

/**
 * This Kernel inserts a block of keys into the hash table.
 * For each keys, it computes a hash position and locates entries
 * with same key within a range of 128 byte after the hash position
 * in parallel.
 *
 * If no match is found, the kernel tries to create a new entry.
 * Therefore, it linear probes the scanned range of 128 byte to find
 * an empty position. If no position has been found, a new hash position
 * will be determined. After MAXTRIALS unsuccessful trials, it inserts
 * the keys into the "noSuccessArea" and continues with the next keys.
 */
template<uint32_t intsPerKey,		// number of ints per key
		uint32_t blocksize,			// number of threads in block
		uint32_t keysPerBlock		// number of keys from bundle processed by this block
>
__global__ void addKernel(
		volatile KeyValuePair<intsPerKey>* table,	// pointer to hash table
		const uint64_t numEntries,			// number of entries in hash table
		const Key<intsPerKey>* keyBundle,	// pointer to keys which are to be inserted
		const uint64_t numKeys,				// number of keys in key bundle
		Key<intsPerKey>* noSuccessArea,		// pointer to area where keys are copied to when they
											// cannot be inserted after MAXTRIALS trials
		uint64_t* numNoSuccessPtr,			// number of unsuccessfully inserted keys
		const uint64_t maxNumNoSuccess		// maximal number of keys in noSuccessArea
		) {

	/**************************************************************************
	 * Declare Variables
	 *************************************************************************/

	// sizes of structures, etc.
	const uint32_t intsPerEntry = intsPerKey + 1;
	const uint32_t entriesPerBlock = blocksize / intsPerEntry;
	const uint32_t counterPosition = intsPerKey;

	// integer-wise interpretation of hash table and kmerBundle
	volatile uint32_t* tableInts = (volatile uint32_t*) table;
	const uint32_t* keyBundleInts = (const uint32_t*) keyBundle;

	// Thread-Indices, etc.
	const uint32_t blockID = blockIdx.x;
	const uint32_t threadID = threadIdx.x;
	const uint32_t subBlockID = threadID % intsPerEntry;
	const uint32_t subBlock = threadID / intsPerEntry;

	// shared variables and pointers to it
	__shared__ uint32_t sharedData[2 + 2 * entriesPerBlock
			+ intsPerKey * keysPerBlock];
	volatile uint32_t* success = &sharedData[0];
	uint32_t* matchOffset = &sharedData[1];
	uint32_t* matchStatus = &sharedData[2];
	uint32_t* lockStatus = &sharedData[2 + entriesPerBlock];
	uint32_t* keyInts = &sharedData[2 + 2 * entriesPerBlock];

	// Local Variables
	uint32_t x;
	uint32_t basePosition;
	uint32_t numEntriesLoaded, probingOffset, numTrials, i;
	uint32_t currentKey;
	volatile uint32_t* addr;
	Key<intsPerKey>* key;

	/**********************************************************************
	 * Load kmers into shared memory
	 *********************************************************************/
	i = threadID;
	while (i < intsPerKey * keysPerBlock && blockID * keysPerBlock * intsPerKey + i < numKeys * intsPerKey) {
		keyInts[i] = keyBundleInts[intsPerKey * blockID * keysPerBlock + i];
#ifdef DEBUG
		printf("blockID=%i threadID=%i loading kmer %i, i=%i, bundlesize=%i\n",
				blockID, threadID, keyInts[i], i, numKeys);
#endif
		i += blocksize;
	}

	/************************************************************************
	 * Try to add keys into hash table
	 ***********************************************************************/

	// for each key in this bundle
	for (currentKey = 0;
			currentKey < keysPerBlock
					&& blockID * keysPerBlock + currentKey < numKeys;
			currentKey++) {

		// Kmer-View on the data in shared memory
		key = (Key<intsPerKey>*) &keyInts[currentKey * intsPerKey];

		if (threadID == 0)
			*success = false;

		/**************************************************************************
		 * Try to add currentKey to hash table
		 *************************************************************************/
		for(numTrials = 0; numTrials < MAXTRIALS && !(*success); numTrials++) {

			// determine probing position
			basePosition = key->hash(numTrials) % numEntries;

#ifdef DEBUG
			printf(
					"blockID=%i, threadID=%i, currentKeyCounter=%i, key[0]=%i, basePosition=%u\n",
					blockID, threadID, currentKey, *((int*) key), basePosition);
#endif

			// number of entries probed in this trial
			numEntriesLoaded = min(entriesPerBlock, (uint32_t) (numEntries - basePosition));

			// thread-own probing address
			addr = &tableInts[basePosition * intsPerEntry + threadID];

			// init flags
			if (threadID < entriesPerBlock) {
				matchStatus[threadID] = true;
			}

			if (threadID == 0)
				*matchOffset = (uint32_t) -1;

#ifdef DEBUG
			printf(
					"blockID=%i, threadID=%i, currentKeyCounter=%i, success=%i, numTrials=%i, basePosition=%i, numEntriesLoaded=%i\n",
					blockID, threadID, currentKey, *success, numTrials,
					basePosition, numEntriesLoaded);
#endif

			// almost all threads participate
			if (threadID < intsPerEntry * numEntriesLoaded) {

				// Coalesced reading of a stripe of memory from hash table:
				// Each Thread reads one integer starting at basePosition
				x = *addr;

#ifdef DEBUG
				printf(
						"blockID=%i, threadID=%i, subwarp=%i, subwarpID=%i: x=%i, subwarpFlags=%i, trial=%i, currentKey=%i\n",
						blockID, threadID, subBlock, subBlockID, x,
						matchStatus[subBlock], numTrials, currentKey);
#endif

				/**************************************************************
				 * See if any position matches currentKey
				 *************************************************************/

				// parallel int-wise comparision with shared kmer
				if (subBlockID < counterPosition) {
					if (x != keyInts[currentKey * intsPerKey + subBlockID])
						matchStatus[subBlock] = false;
				}
				// the last thread of each subwarp checks the counter
				else {

					// determine lock status
					if (x == 0)
						lockStatus[subBlock] = STATUS_FREE;
					else if (x == 0xffffffff)
						lockStatus[subBlock] = STATUS_LOCKED;
					else
						lockStatus[subBlock] = STATUS_OCCUPIED;
				}

				// if match has been found
				if (subBlockID
						== 0 && matchStatus[subBlock] && lockStatus[subBlock] == STATUS_OCCUPIED) {
					*matchOffset = subBlock;
				}


#ifdef DEBUG
				printf(
						"blockID=%i, threadID=%i, subwarp=%i, subwarpID=%i: matchOffset=%i\n",
						blockID, threadID, subBlock, subBlockID, *matchOffset);
#endif

				/**************************************************************
				 * If matching kmer has been found: increase counter
				 *************************************************************/

				// found existing entry in hash table
				if (*matchOffset != (uint32_t) -1) {

					// increase counter
					if (subBlock == *matchOffset
							&& subBlockID == counterPosition) {

						atomicAdd((uint32_t*) addr, 1);
						//*addr = *addr + 1;
#ifdef DEBUG
						printf(
								"blockID = %i: found match at offset = %i (addr=%p)\n",
								blockID, *matchOffset, addr);
#endif

						*success = true;
					}
				}
#ifdef DEBUG
				else {
					if (threadID == 0)
					printf("blockID = %i: no match found\n", blockID);
				}
#endif

				/**************************************************************
				 * If no matching position has been found:
				 * Try to add entry to hash table by linear probing of the
				 * scanned part.
				 **************************************************************/
				for(probingOffset = 0; probingOffset < numEntriesLoaded && !(*success); probingOffset++) {

#ifdef DEBUG
					if (threadID == 0) {
						printf("blockID=%i: probing position %i\n", blockID,
								basePosition + probingOffset);
					}
#endif

					//__threadfence();

					// wait until entry becomes unlocked
					while (lockStatus[probingOffset] == STATUS_LOCKED) {
						
						// re-load table entries from global memory
						x = *addr;

						// init flags
						if (subBlockID == counterPosition)
							matchStatus[subBlock] = true;

#ifdef DEBUG
						printf(
								"blockID=%i, threadID=%i, subwarp=%i, subwarpID=%i: re-loading table address: %p, x=%i, subwarpFlags=%i, trial=%i, currentKey=%i\n",
								blockID, threadID, subBlock, subBlockID, addr, x,
								matchStatus[subBlock], numTrials, currentKey);
#endif

						// parallel int-wise comparision with shared kmer
						if (subBlockID < counterPosition) {
							if (x != keyInts[currentKey * intsPerKey
											+ subBlockID])
								matchStatus[subBlock] = false;
						}
						// the last thread of each subwarp checks the counter
						else {
							// determine lock status
							if (x == 0)
								lockStatus[subBlock] = STATUS_FREE;
							else if (x == 0xffffffff)
								lockStatus[subBlock] = STATUS_LOCKED;
							else
								lockStatus[subBlock] = STATUS_OCCUPIED;
						}
					}

					// entry is now free or occupied

					// entry is occupied by...
					if (lockStatus[probingOffset] == STATUS_OCCUPIED) {

						// ...matching key
						if(matchStatus[probingOffset]) {

							// increase counter by 1
							if (subBlock == probingOffset && subBlockID == counterPosition) {

								atomicAdd((uint32_t*) addr, 1);

#ifdef DEBUG
                	            printf(
           	                 		"blockID = %i: found match at position = %i while reprobing\n",
                                     blockID, basePosition + probingOffset);
#endif

								*success = true;
 							}
						}
						// ...non-matching key
						else {
							// probe next position
						}
					}

					// if entry is free
					else if (lockStatus[probingOffset] == STATUS_FREE) {

						//  try to lock entry
						if (subBlock == probingOffset
								&& subBlockID == counterPosition) {

							if (atomicCAS((uint32_t*) addr, 0, 0xffffffff) == 0) {
								lockStatus[probingOffset] = STATUS_OCCUPIED;
							} else {
								lockStatus[probingOffset] = STATUS_LOCKED;
							}
#ifdef DEBUG
							// if lock was successful
							if (lockStatus[probingOffset] == STATUS_OCCUPIED) {
								printf(
										"blockID=%i, threadID=%i, subwarp=%i, subwarpID=%i, currentKeyCounter=%i: position %i (address %p) was locked\n",
										blockID, threadID, subBlock, subBlockID,
										currentKey,
										basePosition + probingOffset, addr);
							} else {
								printf(
										"blockID=%i, threadID=%i, subwarp=%i, subwarpID=%i, currentKeyCounter=%i: position %i (address %p) was NOT locked\n",
										blockID, threadID, subBlock, subBlockID,
										currentKey,
										basePosition + probingOffset, addr);
							}
#endif
						}

						// if entry could be locked
						if (lockStatus[probingOffset] == STATUS_OCCUPIED) {

							if (subBlock == probingOffset) {

								// copy shared kmer to global memory
								if (subBlockID < counterPosition) {
									*addr = keyInts[currentKey * intsPerKey
											+ subBlockID];
#ifdef DEBUG
									printf(
											"blockID=%i, threadID=%i, subwarp=%i, subwarpID=%i, currentKeyCounter=%i: add entry at position %i (address %p): x=%i\n",
											blockID, threadID, subBlock,
											subBlockID, currentKey,
											basePosition + probingOffset, addr,
											keyInts[currentKey * intsPerKey
											+ subBlockID]);
#endif
								}

								__threadfence();

								if(subBlockID == counterPosition) {

									// set counter to one
									*addr = 1;

									*success = true;
#ifdef DEBUG
									printf(
											"blockID=%i, threadID=%i, subwarp=%i, subwarpID=%i currentKeyCounter=%i, : set counter at position %i (address %p) to 1 -> %i\n",
											blockID, threadID, subBlock,
											subBlockID, currentKey,
											basePosition + probingOffset, addr, *addr);
#endif
								}
							}
						}
						// if entry could not be locked
						else {
							// probe the same entry again
							probingOffset--;
						}
					}
					// error
					else {
						printf("ERROR CASE: %i\n", lockStatus[probingOffset]);
					}
				}

			}
		}
		// if kmer could not be insered after MAXTRIALS trials
		if (!*success) {

			// write kmer to area of shame
			if (threadID == 0) {
#ifdef DEBUG
				printf(
						"blockID=%i, threadID=%i, subwarp=%i, subwarpID=%i:, currentKeyCounter=%i no success after %i trials\n",
						blockID, threadID, subBlock, subBlockID, currentKey,
						numTrials);
#endif

				// get address where to write kmer
				uint64_t numNoSucces = atomicAdd((unsigned long long int*) numNoSuccessPtr, (unsigned long long int) 1);
				if (numNoSucces < maxNumNoSuccess) {
					addr = (uint32_t*) (noSuccessArea + numNoSucces);
					// write kmer to area
					for (uint32_t i = 0; i < intsPerKey; i++)
						addr[i] = keyInts[currentKey * intsPerKey + i];
				}
			}
		}
	}
}

template<uint32_t intsPerKey>
void addToTable(
		KeyValuePair<intsPerKey>* table,	// pointer to the hash table
		const uint64_t numEntries,			// maximum number of entries in table
		const Key<intsPerKey>* keyBundle,	// pointer to the keys that shall be insert
		const uint64_t numKeys,				// number of keys to be inserted
		cudaStream_t& stream,				// cuda stream to use
		Key<intsPerKey>* noSuccessArea,		// pointer to no success area (aka area of shame)
		uint64_t* numNoSuccessPtr,			// pointer to number of unsuccessful inserted keys
		const uint64_t maxNumNoSuccess		// maximum number of keys that can be inserted unsuccessfully
		) {

	// kernel launch parameters
	// number of threads per block is smallest multiple of 32 that is larger than (intsPerKey+1)
	const uint32_t blocksize = ((intsPerKey + 1 + 31) / 32) * 32;
	//const uint32_t processKeysPerBlock = 128 / sizeof(Key<intsPerKey>);
	const uint32_t processKeysPerBlock = 32;
	const uint64_t gridsize = (numKeys + processKeysPerBlock - 1)
			/ processKeysPerBlock;

	/*printf("processKeysPerBlock = %i\n", processKeysPerBlock);
	printf("gridsize            = %lu\n", gridsize);
	printf("numEntries          = %lu\n", numEntries);
	printf("numKeys             = %lu\n", numKeys);
	printf("maxNumNoSuccess     = %lu\n", maxNumNoSuccess);
	printf("sharedMemory        = %lu\n", (sizeof(uint32_t)* (2 + 2 * blocksize / (intsPerKey+1)
	                                      + intsPerKey * processKeysPerBlock)));*/

	// launch kernel
	addKernel
	<intsPerKey, blocksize, processKeysPerBlock>	// template arguments
	<<<gridsize, blocksize, 0, stream>>>			// cuda arguments
	(table, numEntries, keyBundle, numKeys,			// parameters
			noSuccessArea, numNoSuccessPtr, maxNumNoSuccess);

	// check error message
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) {
		throw std::runtime_error(
				std::string(
						std::string("cuda_ds::addToTable::exception1: ")
								+ cudaGetErrorName(err)) + std::string(": ")
						+ std::string(cudaGetErrorString(err)));
	}
}

template<uint32_t intsPerKey>
void sortKeys(Key<intsPerKey>* keys, const uint64_t numKeys) {
	thrust::device_ptr<Key<intsPerKey> > ptr(keys);
	thrust::sort(ptr, ptr + numKeys);
}


/**
 * Export Templates for meaningful ints per key.
 */
#define EXPORT(intsPerKey) 																\
	template void addToTable<intsPerKey>(KeyValuePair<intsPerKey>*, const uint64_t,				\
		const Key<intsPerKey>*, const uint64_t, cudaStream_t&,							\
		Key<intsPerKey>*, uint64_t*, const uint64_t); 									\
	template void sortKeys<intsPerKey>(Key<intsPerKey>*, const uint64_t);

EXPORT(1)
EXPORT(2)
EXPORT(4)
EXPORT(6)
EXPORT(8)
EXPORT(10)
EXPORT(12)
EXPORT(14)
EXPORT(16)
EXPORT(18)
EXPORT(20)
EXPORT(22)
EXPORT(24)
EXPORT(26)
EXPORT(28)
EXPORT(30)
EXPORT(32)
EXPORT(34)

}
}
