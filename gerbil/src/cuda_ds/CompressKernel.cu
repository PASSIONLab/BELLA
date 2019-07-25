#include "../../include/cuda_ds/CountingHashTable/_CountingHashTable.h"

#include <iostream>
#include <stdexcept>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>

namespace cuda_ds {
namespace internal {

/**
 * Determines the first free position i in the hash table with i >= l.
 */
template<uint32_t intsPerKey, uint32_t blocksize>
__device__  inline uint64_t firstFreePosition(KeyValuePair<intsPerKey>* table, const uint64_t numEntries,
		const uint64_t l, const uint64_t r, uint32_t* shared) {

	const uint32_t intsPerEntry = intsPerKey + 1;
	const uint32_t entriesPerBlock = blocksize / intsPerEntry;
	const uint32_t tId = threadIdx.x;

	// Shared Memory and views on it.
	// This is the same physical memory as defined in compressKernel!
	uint32_t* leftPart = &shared[0];
	volatile uint64_t *result = (volatile uint64_t*) &shared[2 * blocksize];

	if (tId == 0)
		*result = (uint64_t) -1;
	__syncthreads();

	uint64_t i = l;
	while (i < r) {

		uint64_t j = i % entriesPerBlock;

		if(j == 0) {

			// read the next stripe of hash table from global memory
			const uint32_t* startAddr = (uint32_t*) (table + i);
			const uint32_t* stopAddr = (uint32_t*) min((uint64_t) (table + i + entriesPerBlock), (uint64_t) (table + numEntries));
			const uint32_t* addr = startAddr + tId;
			if(addr < stopAddr) {
				leftPart[tId] = *addr;
				//printf("thread %u: i=%u loading memory: %u\n", tId, i, leftPart[tId]);
			}
			__syncthreads();
		}

		// linear scan of the part of the hash table
		for (; j < entriesPerBlock && i < r && (*result) == (uint64_t) -1; j++, i++) {
			// only one thread participates
			if (tId == 0) {
				// the first free position is returned
				if (leftPart[j * intsPerEntry + intsPerKey] == 0)
					*result = i;
			}
			__syncthreads();
		}

		if (*result != (uint64_t) -1)
			return *result;
	}

	// no free position found
	return r;
}

/**
 * Determines the last free position i in the hash table with i <= r.
 */
template<uint32_t intsPerKey, uint32_t blocksize>
__device__  inline uint64_t lastOccupiedPosition(KeyValuePair<intsPerKey>* table, const uint64_t numEntries, const uint64_t l,
		const uint64_t r, uint32_t* shared) {

	const uint32_t intsPerEntry = intsPerKey + 1;
	const uint32_t entriesPerBlock = blocksize / intsPerEntry;
	const uint32_t tId = threadIdx.x;

	// Views on the shared memory
	volatile uint32_t* rightPart = &shared[blocksize];
	volatile uint64_t *result = (volatile uint64_t*) &shared[2 * blocksize + 3];

	if (tId == 0)
		*result = (uint64_t) -1;
	__syncthreads();

	// the position from the end of the table
	uint64_t i = numEntries - r - 1;

	// Run to front
	while (numEntries-i-1 > l) {

		uint64_t j = i % entriesPerBlock;

		if(j == 0) {

			const uint64_t numEntriesLoaded = min((uint64_t) entriesPerBlock, numEntries-i);
			const uint32_t* startAddr = (uint32_t*) (table + numEntries - i - numEntriesLoaded);
			const uint32_t* stopAddr = (uint32_t*) (table + numEntries - i);
			const uint32_t* addr = (uint32_t*) (table + numEntries - i - entriesPerBlock) + tId;

			// read the next stripe of hash table from global memory
			if(addr >= startAddr && addr < stopAddr) {
				rightPart[tId] = *addr;
				//printf("thread %u: i=%u loading memory: %u\n", tId, i, rightPart[tId]);
			}
			__syncthreads();
		}

		// linear scan of the part of the hash table
		for (; j < entriesPerBlock && numEntries-i-1 > l && *result == (uint64_t) -1; j++, i++) {
			// only one thread participates
			if (tId == 0) {
				// the last occupied position is returned
				if (rightPart[(entriesPerBlock - j - 1) * intsPerEntry + intsPerKey] != 0)
					*result = numEntries - i -1;
			}
			__syncthreads();
		}

		if (*result != (uint64_t) -1) {
			//printf("r is %u\n", *result);
			return *result;
		}
	}

	// no free position found
	return l;
}

template<uint32_t intsPerKey, uint32_t blocksize>
__device__ inline void my_swap(KeyValuePair<intsPerKey>* table, const uint64_t numEntries,
		const uint64_t l, const uint64_t r, uint32_t* shared) {

	const uint32_t tId = threadIdx.x;
	const uint32_t intsPerEntry = intsPerKey + 1;
	const uint32_t entriesPerBlock = blocksize / intsPerEntry;
	uint32_t* tableInts = (uint32_t*) table;

	//if(tId == 0)
	//	printf("swap %u %u\n", l, r);

	// Views on the shared memory
	uint32_t* leftPart = &shared[0];
	uint32_t* rightPart = &shared[blocksize];

	if (tId < intsPerEntry) {

		// swap in global memory
		const uint64_t x = numEntries - r - 1;
		const uint64_t pos = entriesPerBlock - (x % entriesPerBlock) - 1;
		tableInts[intsPerEntry * l + tId] = rightPart[pos * intsPerEntry + tId];
		tableInts[intsPerEntry * r + tId] = leftPart[(l % entriesPerBlock) * intsPerEntry + tId];
	}
	__syncthreads();

}

/**
 * Compresses the table content to have nonzero entries at the beginning of the table.
 * The number of nonzero entries is stored in *res.
 */
template<uint32_t intsPerKey, uint32_t blocksize>
__global__ void compressKernel(KeyValuePair<intsPerKey>* table,
		const uint64_t numEntries, uint64_t* res) {

	// prepare shared memory to hold parts of the hash table containing l and r
	__shared__ uint32_t shared[2 * blocksize + 4];

	// init left and right marker
	uint64_t l = firstFreePosition<intsPerKey, blocksize>(table, numEntries, 0, numEntries - 1, shared);	// l points to first zero position
	uint64_t r = lastOccupiedPosition<intsPerKey, blocksize>(table, numEntries, 0, numEntries - 1, shared);	// r points to last nonzero position

	/* invariant:
	 *   l points to first zero position.
	 *   r points to last nonzero position.
	 */

	while (l < r) {

		printf("l=%u, r=%u\n", l, r);

		// swap table[l] with table[r])
		my_swap<intsPerKey, blocksize>(table, numEntries, l, r, shared);

		// repair invariant
		l = firstFreePosition<intsPerKey, blocksize>(table, numEntries, l+1, r, shared);
		r = lastOccupiedPosition<intsPerKey, blocksize>(table, numEntries, l, r-1, shared);
	}

	printf("l=%u, r=%u\n", l, r);

	*res = l;
}

template<uint32_t intsPerKey>
uint64_t _compress(KeyValuePair<intsPerKey>* table, const uint64_t numEntries,
		uint64_t* result) {

	const uint32_t blocksize = 1024;

	// launch kernel
	compressKernel<intsPerKey, blocksize> <<<1, blocksize>>>(table, numEntries, result);
	auto err = cudaGetLastError();
	if(err != cudaSuccess) {
		throw std::runtime_error(
				std::string(
						std::string("cuda_ds::compress::exception1: ")
						+ cudaGetErrorName(err)) + std::string(": ")
				+ std::string(cudaGetErrorString(err)));
	}

	uint64_t res;
	cudaMemcpy(&res, result, sizeof(uint64_t), cudaMemcpyDeviceToHost);

	return res;
}

/**
 * Export Templates for meaningful ints per key.
 */
#define EXPORT(x) 																\
	template uint64_t _compress<x>(KeyValuePair<x>*, const uint64_t, uint64_t*);

EXPORT(1)
EXPORT(2)
EXPORT(3)
EXPORT(4)
EXPORT(5)
EXPORT(6)

}
}
