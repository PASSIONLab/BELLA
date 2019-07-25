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

#ifndef KEYVALUEPAIR_HPP_
#define KEYVALUEPAIR_HPP_

#include <string>
#include <sstream>

#include "_Key.h"

namespace cuda_ds {
namespace internal {

/**
 * This Structure encapsulates a key together with a counter.
 */
template<uint32_t intsPerKey>
class KeyValuePair {
public:
	Key<intsPerKey> key;
	uint32_t count;

public:
	__device__ __host__
	KeyValuePair() :
			count(0) {

	}

	/*KeyValuePair(const KeyType& key) :
			key(key), count(0) {

	}*/

	KeyValuePair(const Key<intsPerKey>& key, const uint32_t count) :
			key(key), count(count) {

	}

	const Key<intsPerKey>& getKey() const {
		return key;
	}

	uint32_t getCount() const {
		return count;
	}

	__device__ __host__
	bool operator<(const KeyValuePair<intsPerKey>& kv) const {
		return key < kv.key;
	}

};

}
}

#endif /* KEYVALUEPAIR_HPP_ */
