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

#ifndef KEY_H_
#define KEY_H_

#include <cstdint>

namespace cuda_ds {
namespace internal {

#define C_0  453569
#define C_1  5696063


		/**
 * A generic key class that size can be controlled by the template parameter.
 * @param intsPerKey The number of ints per key.
 */
template<uint32_t intsPerKey>
class Key {

private:
	uint32_t data[intsPerKey];

public:

	/**
	 * Hash function based on byte-wise interpretation of the data
	 */
	__device__ uint64_t hash(const uint32_t trial) const {

		// interpret data as uint16_t
		uint16_t* bytes = (uint16_t*) &data;

		uint64_t h1 = 0;
		for (uint32_t i = 0; i < intsPerKey * 2; i++) {
			h1 += C_0 * h1 + bytes[i];
		}
		return h1 + C_1 * trial * trial;
	}

	/**
	 * Simple int-wise Comparison.
	 */
	__device__ bool operator<(const Key<intsPerKey>& other) const {

		for (uint32_t i = 0; i < intsPerKey; i++) {
			if (data[i] < other.data[i])
				return true;
			else if (data[i] > other.data[i])
				return false;
		}
		return false;
	}

	__device__ __host__ bool operator==(const Key<intsPerKey>& other) const {
		for (uint32_t i = 0; i < intsPerKey; i++) {
			if (data[i] != other.data[i])
				return false;
		}

		return true;
	}

	__device__ __host__ bool operator!=(const Key<intsPerKey>& other) const {
		return !operator ==(other);
	}

};

}
}

#endif /* KEY_H_ */
