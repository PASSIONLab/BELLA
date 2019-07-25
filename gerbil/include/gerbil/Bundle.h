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

#ifndef BUNDLE_H_
#define BUNDLE_H_

#include "KMer.h"

namespace gerbil {


//#define USE_MEM_COPY

	class FastBundle {
	public:
		char data[FAST_BUNDLE_DATA_SIZE_B + 1];
		uint32 size;
		TFileCompr compressType;

		FastBundle() :
				size(0) {
		}

		bool isFull() {
			return size + FAST_BLOCK_SIZE_B > FAST_BUNDLE_DATA_SIZE_B;
		}

		void inline finalize(const TFileCompr& pFileCompr) {
			data[size] = '\0';
			compressType = pFileCompr;
		}

		void inline clear() {
			size = 0;
		}
	};

/*
 * Bundle for single reads
 */
	class ReadBundle {
		static uint K;
	public:
		byte *data;                // data (header + reads)
		uint32 *readsCount;        // number of reads
		uint32 *readOffsets;// offset of reads in data, mark start (inclusive) and end (exclusive)

		ReadBundle();

		~ReadBundle();

		bool isEmpty() const;

		// reset ReadBundle
		void clear();

		// adds a single read to a ReadBundle
		bool add(const uint32 &length, char *read);

		// expands last read, returns false on failure
		bool expand(const uint32 &expLength, char *expRead);

		// transfers last read from this to readBundle
		// Cond.: readsCount > 0
		bool transfer(ReadBundle *readBundle);

		// transfers the last k-1 bases from the last of this bundle to readBundle
		// Cond.: readsCount > 0
		bool transferKm1(ReadBundle *readbundle);

		void print();

		static void setK(uint k) { K = k; }
	};

/*
 * Bundle for s-mers
 */
	class SuperBundle {

	private:
		byte *_next;
		bool _finalized;
	public:
		byte data[SUPER_BUNDLE_DATA_SIZE_B];                // data
		uint_tfn tempFileId;                                // id of tempFile
		uint_tfn tempFileRun;                               // run of TempFile

		uint64 sMerNumber;                                    // number of s-mers
		uint64 kMerNumber;                                    // number of k-mers

		SuperBundle();

		~SuperBundle();

		void finalize();

		void clear();

		bool add(const byte *smer, const uint16 &length, const uint32_t &k);

		bool next(byte *&smer, uint16 &length);

		bool isEmpty() const;

		bool merge(const SuperBundle &sb);

//debugging, working only with "add"
		uint32 getSize() {
			return _next - data + 1;
		}
	};

/*
 * bundle for k-mers
 */
	template<uint32_t K, uint32_t bufferSize>
	class KMerBundle {
	private:
		KMer<K> *_data;            // data (k-mers)
		KMer<K> *_next;            // next place for reading k-mers
		KMer<K> *_last;            // next place for writing k-mers
		KMer<K> *_end;            // end of data (exclusive)
		uint_tfn _tempFileId;    // id of TempFile
		uint_tfn _tempFileRun;  // run of TempFile

		static constexpr size_t dataSize() { return bufferSize / sizeof(KMer<K>); }

		static constexpr size_t dataSize_B() { return dataSize() * sizeof(KMer<K>); }

	public:
		KMerBundle() :
				_tempFileId(TEMPFILEID_NONE), _tempFileRun(0) {
			size_t l = bufferSize / sizeof(KMer<K>);
			_data = new KMer<K>[dataSize()];
			_last = _next = _data;
			_end = _data + dataSize();
		}

		~KMerBundle() {
			delete[] _data;
		}

		inline bool isEmpty() const {
			return _last == _data;
		}

		inline void setTempFileId(const uint_tfn &binId) {
			_tempFileId = binId;
		}

		inline void setTempFileRun(const uint_tfn &run) {
			_tempFileRun = run;
		}

		inline uint_tfn getTempFileId() const {
			return _tempFileId;
		}

		inline uint_tfn getTempFileRun() const {
			return _tempFileRun;
		}

		inline KMer<K> *getData() {
			return _data;
		}

		inline void clear() {
			_next = _data;
			_last = _data;
			_tempFileId = TEMPFILEID_NONE;
			_tempFileRun = 0;
		}

		inline void store(FILE *&file) {
			assert(_last == _end);
			fwrite((char *) _data, 1, dataSize_B(), file);
		}

		inline void load(FILE *&file) {
			clear();
			if (fread((char *) _data, 1, dataSize_B(), file) != dataSize_B()) {
				std::cerr << "ERROR: reload of k-mers failed\n";
				exit(7);
			}
			_last = _end;
		}

		inline bool add(const KMer<K> &kMer) {
			if (_last < _end) {
				(_last++)->set(kMer);
				return true;
			}
			return false;
		}

		inline uint32_t count() const {
			return (uint32_t) (((uint64_t) _last - (uint64_t) _data) / sizeof(KMer<K>));
		}

		inline bool next(KMer<K> *&kMer) {
			if (_next < _last) {
				kMer = (_next++);
				return true;
			}
			return false;
		}

		inline void copyAndInc(KMer<K> *&kMer) {
			while (_next < _last)
				(kMer++)->set(*(_next++));
		}

		void print() const {
			const KMer<K> *i = _data;
			while (i < _last)
				printKMer(*(i++));
		}
	};

// declare KmerBundles for CPU
	namespace cpu {
		template<uint32_t K>
		class KMerBundle : public ::gerbil::KMerBundle<K, KMER_BUNDLE_DATA_SIZE_B> {
		};
	}

// declare KmerBundles for GPU
	namespace gpu {
		template<uint32_t K>
		class KMerBundle : public ::gerbil::KMerBundle<K, GPU_KMER_BUNDLE_DATA_SIZE> {
		};
	}

/*
 * bundle for k-mer counts
 */
	class KmcBundle {
		byte *_data;            // data
		byte *_next;            // next free place
	public:
		KmcBundle();

		~KmcBundle();

		// reset KmcBundle
		void clear();

		// add k-mer with count
		template<unsigned K>
		bool add(const KMer<K> &kMer, const uint32 &val);

		bool isEmpty() const;

		uint32 getSize() const;

		const byte *getData() const;
	};

	template<unsigned K>
	bool KmcBundle::add(const KMer<K> &kMer, const uint32 &val) {
		if (_next + sizeof(KMer<K>) + 1
		    + (val < 255 ? 1 : 5) >= _data + KMC_BUNDLE_DATA_SIZE_B)
			return false;
		if (val < 255)    // 1 byte for counts < 255
			*(_next++) = val;
		else {            // 1 byte (Count greater than 254) + 4 bytes for count
			*(_next++) = 0xff;
			*((uint32 *) _next) = val;
			_next += 4;
		}
		// add kMer after count
		kMer.toByte(_next);
		_next += getKMerCompactByteNumbers<K>();
		return true;
	}

	inline bool KmcBundle::isEmpty() const {
		return _next == _data;
	}

	inline uint32 KmcBundle::getSize() const {
		return _next - _data;
	}

	inline const byte *KmcBundle::getData() const {
		return _data;
	}

}

#endif /* BUNDLE_H_ */
