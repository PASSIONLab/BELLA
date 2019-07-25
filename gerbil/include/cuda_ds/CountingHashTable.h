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

#ifndef COUNTINGHASHTABLE_H_
#define COUNTINGHASHTABLE_H_

#include "CountingHashTable/_CountingHashTable.h"

namespace cuda_ds {

//#define internal::intsPerKey<KeyType>() sizeof(KeyType) / sizeof(uint32_t)

/**
 * A Data Structure for Counting keys. It is based on the idea
 * of a counting hash table.
 * To add insert inserts a bundle of keys into the hash table.
 * The key counts can be determined by calling the count method.
 *
 * @param KeyType The type of the keys. Attention! KeyType's size
 * has to be a multiple of sizeof(int) and has to be smaller than
 * 32*sizeof(int).
 * @param keyBufferSize Insert the keys in chunks of size keyBufferSize.
 */
template<class KeyType, uint64_t keyBufferSize = 512 * 1024>
class CountingHashTable: public internal::CountingHashTable<
		internal::intsPerKey<KeyType>(), keyBufferSize> {

public:

	/**
	 * Constructor. Just a Wrapper to the internal type.
	 */
	CountingHashTable(const uint32_t devID = 0) :
			internal::CountingHashTable<internal::intsPerKey<KeyType>(),
					keyBufferSize>(devID) {

		// check whether size of KeyType is a multiple of sizeof(uint32_t)
		if (sizeof(KeyType) % sizeof(uint32_t) != 0) {
			std::runtime_error(
					std::string(
							"cuda_ds::CountingHashTable::exception: size of keytype is not a multiple of sizeof(int)."));
		}
	}

	/**
	 * Insert a bundle of keys to the hash table.
	 * @param keyBundle A bundle of keys.
	 * @param numKeys The number of keys to insert into the table.
	 * @return True, if all keys could be insert safely.
	 *         False, if there is not enough room to insert this
	 *         number of keys safely. In this case, no keys are
	 *         insert. An extract is required.
	 */
	bool insert(const KeyType* keyBundle, const uint64_t numKeys) {

		// cast to internal Key type and call internal add
		return this->add(
				(const internal::Key<internal::intsPerKey<KeyType>()>*) keyBundle,
				numKeys);
	}

	/**
	 * Wait until all keys are inserted and count the result.
	 */
	template<uint64_t countBufferSize = 16 * 1024 * 1024>
	void extract(std::vector<std::pair<KeyType, uint32_t>>& counts) const {

		// set device id
		cudaSetDevice(this->devID);

		// wait for threads to success
		/*cudaError_t err;
		 err = cudaDeviceSynchronize();
		 if (err != cudaSuccess) {
		 throw std::runtime_error(
		 std::string(
		 std::string("cuda_ds::exception: ")
		 + cudaGetErrorName(err)) + std::string(": ")
		 + std::string(cudaGetErrorString(err)));
		 }*/

		// compress table. this is a bad idea (worse performance)
		//uint32_t numDifferentKeys = this->compressEntries();
		//printf("keys in hash table: %u\n", numDifferentKeys);
		// Load table to host memory in small chunks
		uint64_t itemsPerBuffer =
				countBufferSize
						/ sizeof(internal::KeyValuePair<
								internal::intsPerKey<KeyType>()>);
		internal::KeyValuePair<internal::intsPerKey<KeyType>()> *hostTable =
				new internal::KeyValuePair<internal::intsPerKey<KeyType>()>[itemsPerBuffer];

		for (uint32_t i = 0; i * itemsPerBuffer < this->getNumEntries(); i++) {

			uint64_t itemsToFetch =
					(i + 1) * itemsPerBuffer < this->getNumEntries() ?
							itemsPerBuffer :
							this->getNumEntries() % itemsPerBuffer;

			cudaMemcpy(hostTable, this->table + i * itemsPerBuffer,
					itemsToFetch
							* sizeof(internal::KeyValuePair<
									internal::intsPerKey<KeyType>()>),
					cudaMemcpyDeviceToHost);

			// copy kmer counts to result vector
			for (int j = 0; j < itemsToFetch; j++) {
				if (hostTable[j].getCount() > 0) {
					//std::pair<KeyType, uint32_t> kv = std::make_pair(*((KeyType*) &hostTable[j].getKey()), hostTable[j].getCount());
					//counts.push_back(kv);
					counts.push_back(
							*((std::pair<KeyType, uint32_t>*) &hostTable[j]));
				}
			}
		}

		delete[] hostTable;

		// Determine if there are any unsuccessfully inserted kmers
		uint64_t numNoSuccess;
		cudaMemcpy(&numNoSuccess, this->numNoSuccessPtr, sizeof(uint64_t),
				cudaMemcpyDeviceToHost);

		printf("keys in free area: %u\n", numNoSuccess);

		assert(numNoSuccess <= this->maxNumNoSuccess);

		// if there are some
		if (numNoSuccess > 0) {

			// sort area of shame
			sortKeys(this->noSuccessArea, numNoSuccess);

			// Load keys (in small chunks)
			internal::Key<internal::intsPerKey<KeyType>()> *hostArea;
			itemsPerBuffer = countBufferSize
					/ sizeof(internal::Key<internal::intsPerKey<KeyType>()>);
			hostArea =
					new internal::Key<internal::intsPerKey<KeyType>()>[itemsPerBuffer];

			// load first buffer
			uint32_t itemsToFetch =
					itemsPerBuffer <= numNoSuccess ?
							itemsPerBuffer : numNoSuccess % itemsPerBuffer;

			cudaMemcpy(hostArea, this->noSuccessArea,
					itemsToFetch
							* sizeof(internal::Key<
									internal::intsPerKey<KeyType>()>),
					cudaMemcpyDeviceToHost);

			// compare with last seen key
			internal::Key<internal::intsPerKey<KeyType>()> lastKey = hostArea[0];
			uint32_t count = 1;

			// compress and add to result vector
			for (uint32_t j = 1; j < itemsToFetch; j++) {
				if (hostArea[j] != lastKey) {
					counts.push_back(
							std::make_pair(*(KeyType*) &lastKey, count));
					lastKey = hostArea[j];
					count = 1;
				} else
					count++;
			}

			// load other buffers (if any)
			for (uint32_t i = 1; i * itemsPerBuffer < numNoSuccess; i++) {

				itemsToFetch =
						(i + 1) * itemsPerBuffer < numNoSuccess ?
								itemsPerBuffer : numNoSuccess % itemsPerBuffer;

				cudaMemcpy(hostArea, this->noSuccessArea + i * itemsPerBuffer,
						itemsToFetch
								* sizeof(internal::Key<
										internal::intsPerKey<KeyType>()>),
						cudaMemcpyDeviceToHost);

				// compress and add to result vector
				for (uint32_t j = 0; j < itemsToFetch; j++) {
					if (hostArea[j] != lastKey) {
						counts.push_back(
								std::make_pair(*(KeyType*) &lastKey, count));
						lastKey = hostArea[j];
						count = 1;
					} else
						count++;
				}
			}

			// insert last seen kmer
			counts.push_back(std::make_pair(*(KeyType*) &lastKey, count));

			delete[] hostArea;
		}

	}

	/**
	 * Print the table.
	 */
	std::string toString() const {
		std::stringstream ss;
		ss << "Table:\n";
		internal::KeyValuePair<internal::intsPerKey<KeyType>()> *hostTable =
				new internal::KeyValuePair<internal::intsPerKey<KeyType>()>[this->numEntries];

		cudaMemcpy(hostTable, this->table,
				this->numEntries
						* sizeof(internal::KeyValuePair<
								internal::intsPerKey<KeyType>()>),
				cudaMemcpyDeviceToHost);

		for (int i = 0; i < this->numEntries; i++) {
			ss << *((KeyType*) &hostTable[i].getKey()) << ": "
					<< hostTable[i].getCount() << "\n";
		}
		delete[] hostTable;
		return ss.str();
	}

};

}

template<class KeyType, uint64_t keyBufferSize>
std::ostream &operator<<(std::ostream &os,
		cuda_ds::CountingHashTable<KeyType, keyBufferSize> const & table) {
	return os << table.toString();
}

#endif /* COUNTINGHASHTABLE_H_ */
