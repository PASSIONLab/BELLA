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

#ifndef KMER_H_
#define KMER_H_

#include "types.h"
#include <atomic>
#include <xmmintrin.h>
#include <algorithm>

namespace gerbil {

// number of bytes for k-mer
#define GET_KMER_B(k) (k / 4 + 1)

// k-mer type, use of 32 Bit structures if b <= 4, else 64 Bit structures
#define GET_KMER_T(b) (b <= 4 ? 4 : 8)

// number of necessary t-bytes-structures to store a single k-mer
#define GET_KMER_C(b, t) ((b - 1) / t + 1)

// empty k-mer, k > 4  (T == 8)
#define KMER_EMPTY_T8 0x0000000000000001

// empty k-mer, k <= 4 (T == 4)
#define KMER_EMPTY_T4 0x00000001

	/*
	 * empty template
	 * general structure of KMers
	 */
	template<unsigned K, unsigned B = GET_KMER_B(K), unsigned T = GET_KMER_T(B),
			unsigned C = GET_KMER_C(B, T)>
	struct KMer {
	};

	/*
	 * KMer for k >= 32
	 * --> T == 8
	 * --> C >  1
	 */
	template<unsigned K, unsigned B, unsigned C>
	struct KMer<K, B, 8, C> {
		uint64_t data[C];
		static constexpr uint64_t _c_offset = 64 - (2 * (K % 32));
		static constexpr uint64_t _c_mask = 0xffffffffffffffff << _c_offset;
	public:
		/*
		 * clears this KMer
		 */
		inline void clear() {
			*data = KMER_EMPTY_T8;
		}

		/*
		 * copy function
		 */
		inline void set(const KMer<K> &kMer) {
			for (uint c = 0; c < C; ++c)
				data[c] = kMer.data[c];
		}

		/*
		 * fills the data of this KMer in dependence of the passed bytes (== sequence of bases, e.g. s-mer)
		 */
		inline void set(const byte *const &bytes) {
			const byte *cb(bytes);
			for (uint c(C); --c;) {
				uint64 &cd(data[c]);
				cd = *cb;
				cd <<= 8;
				cd |= *(++cb);
				cd <<= 8;
				cd |= *(++cb);
				cd <<= 8;
				cd |= *(++cb);
				cd <<= 8;
				cd |= *(++cb);
				cd <<= 8;
				cd |= *(++cb);
				cd <<= 8;
				cd |= *(++cb);
				cd <<= 8;
				cd |= *((++cb)++);
			}
			// first block (some unused space)
			uint64 &cd(*data);
			cd = *cb;
			cd <<= 8;
			cd |= *(++cb);
			cd <<= 8;
			cd |= *(++cb);
			cd <<= 8;
			cd |= *(++cb);
			cd <<= 8;
			cd |= *(++cb);
			cd <<= 8;
			cd |= *(++cb);
			cd <<= 8;
			cd |= *(++cb);
			cd <<= 8;
			cd |= *(++cb);
			// cleaning up unused space
			cd &= _c_mask;
		}

		/*
		 * fills the data of this KMer in dependence of the passed bytes (== sequence of bases, e.g. s-mer)
		 * INVERSE k-mer
		 */
		inline void setInv(const byte *const &bytes) {
			uint64 *p;
			for (uint c = 0; c < C; ++c) {
				// inverts the bases and swaps every two adjacent 8-bit blocks
				*(p = data + c) = ~((((((((((((((((uint64) bytes[8 * c + 7]) << 8)
				                                | bytes[8 * c + 6]) << 8) | bytes[8 * c + 5]) << 8)
				                            | bytes[8 * c + 4]) << 8) | bytes[8 * c + 3]) << 8)
				                        | bytes[8 * c + 2]) << 8) | bytes[8 * c + 1]) << 8)
				                    | bytes[8 * c + 0]);
				// swaps every two adjacent 4-bit blocks
				*p = ((*p & 0xf0f0f0f0f0f0f0f0) >> 4)
				     | ((*p & 0x0f0f0f0f0f0f0f0f) << 4);
				//swaps every two adjacent 2-bit blocks (== single bases)
				*p = ((*p & 0xcccccccccccccccc) >> 2)
				     | ((*p & 0x3333333333333333) << 2);
			}
			if (K % 32) { // just ignore warnings...
#pragma warning(push, 0)
				for (uint c = C; --c;) {
					data[c] <<= _c_offset;
					data[c] |= (data[c - 1] >> (2 * (K % 32)));
				}
				data[0] <<= _c_offset;
#pragma warning(pop)
			} else {
				for (uint c = C; --c;)
					data[c] = data[c - 1];
				*data = 0;
			}
		}

		/*
		 * sets the k-mer and the inverse k-mer
		 */
		static inline void set(const byte *const &bytes, KMer<K> &kMer,
		                       KMer<K> &iKMer) {
			kMer.set(bytes);
			iKMer.setInv(bytes);
		}

		/*
		 * returns a reference to the normalized KMer
		 */
		inline const KMer<K> &getNormalized(const KMer<K> &invKMer) const {
			for (uint i = C - 1; i; i--)
				if (data[i] < invKMer.data[i])
					return *this;
				else if (data[i] > invKMer.data[i])
					return invKMer;
			return data[0] < invKMer.data[0] ? *this : invKMer;
		}

		/*
		 * shifts the whole k-mer and adds a new base
		 */
		inline void next(const uint8_t &b) {
			if (K % 32) {
				for (uint c = C - 1; c > 0; --c)
					data[c] = (data[c] << 2) | (data[c - 1] >> 62);
				data[0] = (data[0] << 2) | (((uint64) b) << _c_offset);
			} else {
				// entire first block is unused
				for (uint c = C - 1; c > 1; --c)
					data[c] = (data[c] << 2) | (data[c - 1] >> 62);
				data[1] = (data[1] << 2) | b;
			}
		}

		/*
		 * shifts the whole inverse k-mer and adds a new base
		 */
		inline void nextInv(const uint8_t &b) {
			data[0] = ((data[0] >> 2) | ((data[1] & 0x3) << 62)) & _c_mask;
			for (uint c = 1; c < C - 1; ++c)
				data[c] = (data[c] >> 2) | ((data[c + 1] & 0x3) << 62);
			data[C - 1] >>= 2;
			data[C - 1] |= ((uint64) (~b)) << 62;
		}

		/*
		 * checks, whether data are equal or not
		 */
		inline bool isEqual(const KMer<K> &kMer) const {
			if (*data != *kMer.data)
				return false;
			for (uint c = 1; c < C; ++c)
				if (data[c] != kMer.data[c])
					return false;
			return true;
		}

		/*
		 * checks, whether it is empty
		 * DEPRECATED
		 */
		inline bool isEmpty() const {
			return *data == KMER_EMPTY_T8;
		}

		/*
		 * magic hash function for counting
		 */
		inline uint64_t getHash() const {
			/*if (C == 2)
				return ((data[0] + data[1]) * C_0 + C_1);
			else if(C == 3)
				return ((data[0] + data[1] + data[2]) * C_0 + C_1);
			else {
				uint64_t x = (((uint64) data[0] << (K % 32)) ^ data[1]) * C_0 + C_1;
				for (uint c = 2; c < C; ++c)
					x ^= (data[c] << ((K % 32) * (c & 0x1)));
				return x;
			}*/

			uint64_t x = 0;
			for(int c = 0; c<C; ++c) {
				x += C_0 * x + data[c];
			}
			return x;
			//return std::accumulate(data, data+C, (uint64_t) 0) * C_0;
		}

		/*
		 * second magic hash function for distribution on threads
		 */
		inline uint64_t getPartHash() const {


			const uint16_t * bytes = (uint16_t*) (data + (C == 2 ? 0 : 1));
			uint64_t x = 0;
			for(int i=0; i<8; i++) {
				x += C_0 * x + bytes[i];
			}
			return x;

			//return (data[1] >> 1) ^ data[C == 2 ? 0 : 2];
			//return (data[1] * C_1 + data[C == 2 ? 0 : 2]) * C_2;
		}

		/*
		 * opposite of set
		 * returns a simple byte representation of this KMer
		 */
		inline void toByte(byte *const bytes) const {
			byte *bs(bytes);
			for (uint64_t c(C); --c;) {
				const uint64 &cd(data[c]);
				*bs = cd >> 56;
				*(++bs) = cd >> 48;
				*(++bs) = cd >> 40;
				*(++bs) = cd >> 32;
				*(++bs) = cd >> 24;
				*(++bs) = cd >> 16;
				*(++bs) = cd >> 8;
				*((++bs)++) = cd;
			}
			const uint64 &cd(*data);
			if ((K & 0x1f) > 0) {
				*bs = cd >> 56;
				if ((K & 0x1f) > 4) {
					*(++bs) = cd >> 48;
					if ((K & 0x1f) > 8) {
						*(++bs) = cd >> 40;
						if ((K & 0x1f) > 12) {
							*(++bs) = cd >> 32;
							if ((K & 0x1f) > 16) {
								*(++bs) = cd >> 24;
								if ((K & 0x1f) > 20) {
									*(++bs) = cd >> 16;
									if ((K & 0x1f) > 24) {
										*(++bs) = cd >> 8;
										if ((K & 0x1f) > 28) {
											*(++bs) = cd;
										}
									}
								}
							}
						}
					}
				}
			}
		}

		/*
		 * returns true if the data are equal
		 */
		inline bool operator==(const KMer<K, B, 8, C> &other) const {
			for (int i = 0; i < C; i++) {
				if (data[i] != other.data[i])
					return false;
			}
			return true;
		}

		/*
		 * returns true if the data are not equal
		 */
		inline bool operator!=(const KMer<K, B, 8, C> &other) const {
			return !operator==(other);
		}
	};


/*
 * KMer for 16 <= k <= 31
 * --> 5 <= B <= 8
 * --> T == 8
 * --> C == 1
 */
	template<unsigned K, unsigned B>
	struct KMer<K, B, 8, 1> {
		uint64_t data;
		static constexpr uint64_t _c_offset = 64 - (2 * K);
		static constexpr uint64_t _c_mask = 0xffffffffffffffff << _c_offset;
	public:
		/*
		 * clears this KMer
		 */
		inline void clear() {
			data = KMER_EMPTY_T8;
		}

		/*
		 * simple copy function
		 */
		inline void set(const KMer<K> &kMer) {
			data = kMer.data;
		}


		/*
		 * fills the data of this KMer in dependence of the passed bytes (== sequence of bases, e.g. s-mer)
		 */
		inline void set(const byte *const &bytes) {
			data = *bytes;
			data <<= 8;
			data |= bytes[1];
			data <<= 8;
			data |= bytes[2];
			data <<= 8;
			data |= bytes[3];
			if (B >= 5) {
				data <<= 8;
				data |= bytes[4];
				if (B >= 6) {
					data <<= 8;
					data |= bytes[5];
					if (B >= 7) {
						data <<= 8;
						data |= bytes[6];
						data <<= 8;
						if (B == 8)
							data |= bytes[7];
					} else
						data <<= 16;
				} else
					data <<= 24;
			} else
				data <<= 32;
			data &= _c_mask;
		}

		/*
		 * fills the data of this KMer in dependence of the passed bytes (== sequence of bases, e.g. s-mer)
		 * INVERSE k-mer
		 */
		inline void setInv(const KMer<K> &kMer) {
			data = ~kMer.data;
			data = (data >> 32) | (data << 32);
			data = ((data >> 16) & 0x0000ffff0000ffff)
			       | ((data & 0x0000ffff0000ffff) << 16);
			data = ((data >> 8) & 0x00ff00ff00ff00ff)
			       | ((data & 0x00ff00ff00ff00ff) << 8);
			data = ((data >> 4) & 0x0f0f0f0f0f0f0f0f)
			       | ((data & 0x0f0f0f0f0f0f0f0f) << 4);
			data = ((data >> 2) & 0x3333333333333333)
			       | ((data & 0x3333333333333333) << 2);
			data <<= _c_offset;
		}


		/*
		 * sets the k-mer and the inverse k-mer
		 */
		static inline void set(const byte *const &bytes, KMer<K> &kMer,
		                       KMer<K> &iKMer) {
			kMer.set(bytes);
			iKMer.setInv(kMer);
		}

		/*
		 * shifts the whole k-mer and adds a new base
		 */
		inline void next(const uint8_t &b) {
			data <<= 2;
			data |= (static_cast<uint64>(b) << _c_offset);
		}

		/*
		 * shifts the whole inverse k-mer and adds a new base
		 */
		inline void nextInv(const uint8_t &b) {
			data >>= 2;
			data |= (~static_cast<uint64>(b)) << 62;
			data &= _c_mask;
		}

		/*
		 * returns a reference to the normalized KMer
		 */
		inline const KMer<K> &getNormalized(const KMer<K> &invKMer) const {
			return data < invKMer.data ? *this : invKMer;
		}

		/*
		 * checks, whether data are equal or not
		 */
		inline bool isEqual(const KMer<K> &kMer) const {
			return data == kMer.data;
		}

		/*
		 * checks, whether it is empty
		 * DEPRECATED
		 */
		inline bool isEmpty() const {
			return data == KMER_EMPTY_T8;
		}

		/*
		 * magic hash function for counting
		 */
		inline uint64_t getHash() const {
			return ((data >> _c_offset) * C_0 + C_1) % 849399569653;
		}

		/*
		 * second magic hash function for distribution on threads
		 */
		inline uint64_t getPartHash() const {
			return ((data >> _c_offset) * C_3 + C_4) % 849399569653;
			//uint64_t x = (data >> _c_offset) & 0x5555555555555555;
			//return x | (x >> ((K - 1) | 0x1));
		}

		/*
		 * opposite of set
		 * returns a simple byte representation of this KMer
		 */
		inline void toByte(byte *const bytes) const {
			bytes[0] = data >> 56;
			bytes[1] = data >> 48;
			bytes[2] = data >> 40;
			bytes[3] = data >> 32;
			if (B >= 5)
				bytes[4] = data >> 24;
			if (B >= 6)
				bytes[5] = data >> 16;
			if (B >= 7)
				bytes[6] = data >> 8;
			if (B == 8)
				bytes[7] = data;
		}

		/*
		 * returns true if the data are equal
		 */
		inline bool operator==(const KMer<K, B, 8, 1> &other) const {
			return data == other.data;
		}

		/*
		 * returns true if the data are not equal
		 */
		inline bool operator!=(const KMer<K, B, 8, 1> &other) const {
			return !operator==(other);
		}
	};

	/*
	 * KMer for k <= 15
	 * --> B <= 4
	 * --> T == 4
	 * --> C == 1
	 */
	template<unsigned K, unsigned B>
	struct KMer<K, B, 4, 1> {
		uint32 data;
		static constexpr uint32 _c_offset = 32 - (2 * K);
		static constexpr uint32 _c_mask = 0xffffffff << _c_offset;
	public:
		/*
		 * clears this KMer
		 */
		inline void clear() {
			data = KMER_EMPTY_T4;
		}

		/*
		 * simple copy function
		 */
		inline void set(const KMer<K> &kMer) {
			data = kMer.data;
		}

		/*
		 * fills the data of this KMer in dependence of the passed bytes (== sequence of bases, e.g. s-mer)
		 */
		inline void set(const byte *const &bytes) {
			data = *bytes;
			data <<= 8;
			data |= bytes[1];
			data <<= 8;
			data |= bytes[2];
			data <<= 8;
			data |= bytes[3];
			data &= _c_mask;
		}

		/*
		 * fills the data of this KMer in dependence of the passed bytes (== sequence of bases, e.g. s-mer)
		 * INVERSE k-mer
		 */
		inline void setInv(const KMer<K> &kMer) {
			data = ~kMer.data;
			data = (data >> 16) | (data << 16);
			data = ((data >> 8) & 0x00ff00ff) | ((data & 0x00ff00ff) << 8);
			data = ((data >> 4) & 0x0f0f0f0f) | ((data & 0x0f0f0f0f) << 4);
			data = ((data >> 2) & 0x33333333) | ((data & 0x33333333) << 2);
			data <<= _c_offset;
		}

		/*
		 * sets the k-mer and the inverse k-mer
		 */
		static inline void set(const byte *const &bytes, KMer<K> &kMer,
		                       KMer<K> &iKMer) {
			kMer.set(bytes);
			iKMer.setInv(kMer);
		}

		/*
		 * shifts the whole k-mer and adds a new base
		 */
		inline void next(const uint8_t &b) {
			data <<= 2;
			data |= (static_cast<uint32>(b) << _c_offset);
		}

		/*
		 * shifts the whole inverse k-mer and adds a new base
		 */
		inline void nextInv(const uint8_t &b) {
			data >>= 2;
			data |= (~static_cast<uint32>(b)) << 30;
			data &= _c_mask;
		}

		/*
		 * returns a reference to the normalized KMer
		 */
		inline const KMer<K> &getNormalized(const KMer<K> &invKMer) const {
			return data < invKMer.data ? *this : invKMer;
		}

		/*
		 * checks, whether data are equal or not
		 */
		inline bool isEqual(const KMer<K> &kMer) const {
			return data == kMer.data;
		}

		/*
		 * checks, whether it is empty
		 * DEPRECATED
		 */
		inline bool isEmpty() const {
			return data == KMER_EMPTY_T4;
		}

		/*
		 * magic hash function for counting
		 */
		inline uint64_t getHash() const {
			return ((uint64) (data >> _c_offset) * C_0 + C_1) % C_2;
		}

		/*
		 * second magic hash function for distribution on threads
		 */
		inline uint64_t getPartHash() const {
			uint64_t x((data >> _c_offset) & 0x55555555);
			return x | (x >> ((K - 1) | 0x1));
		}

		/*
		 * opposite of set
		 * returns a simple byte representation of this KMer
		 */
		inline void toByte(byte *const bytes) const {
			bytes[0] = data >> 24;
			bytes[1] = data >> 16;
			bytes[2] = data >> 8;
			bytes[3] = data;
		}

		/*
		 * returns true if the data are equal
		 */
		inline bool operator==(const KMer<K, B, 4, 1> &other) const {
			return data == other.data;
		}

		/*
		 * returns true if the data are not equal
		 */
		inline bool operator!=(const KMer<K, B, 4, 1> &other) const {
			return !operator==(other);
		}
	};

	template<unsigned K, unsigned B, unsigned C>
	bool operator<(const KMer<K, B, 8, C> &l, const KMer<K, B, 8, C> &r) {
		for (uint i(1); i < C; ++i)
			if (l.data[i] < r.data[i])
				return true;
		return *(l.data) < *(r.data);
	}

	template<unsigned K, unsigned B>
	bool operator<(const KMer<K, B, 8, 1> &l, const KMer<K, B, 8, 1> &r) {
		return l.data < r.data;
	}

	template<unsigned K, unsigned B>
	bool operator<(const KMer<K, B, 4, 1> &l, const KMer<K, B, 4, 1> &r) {
		return l.data < r.data;
	}

/*
 * returns the minimum number of bytes to save this KMer
 * --> 4 bpB
 */
	template<unsigned K>
	constexpr uint16 getKMerCompactByteNumbers() {
		return (K + 3) >> 2;
	}

/*
 * prints a readable representation of this KMer
 * DEBUG ONLY!
 */
	template<unsigned K>
	void printKMer(const KMer<K> &kMer) {
		byte *a = new byte[sizeof(KMer<K>)];
		kMer.toByte(a);
		printByteCodedSeqN(a, K);
		delete[] a;
	}

/*
 * equal to sizeof(KMer<K>)
 * DEPRECATED
 */
	inline uint32 getKMerByteNumbers(const uint32_t &k) {
		uint32 b = GET_KMER_B(k);
		uint32 t = GET_KMER_T(b);
		return GET_KMER_C(b, t) * t;
	}

/*
 * equal to getKMerCompactByteNumbers<K>
 */
	inline uint32_t getKMerCompactByteNumbers(const uint32_t &k) {
		return (k + 3) / 4;
	}

}

#endif /* KMER_H_ */
