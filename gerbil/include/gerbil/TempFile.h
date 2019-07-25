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

#ifndef BINFILE_H_
#define BINFILE_H_

#include "Bundle.h"

namespace gerbil {

	class TempFile {
		static uint_tfn __nextId;       // guaranteed unique ids

		uint_tfn _id;                   // id of bin file

		std::string _filename;          // filename
		FILE *_file;                    // file
		uint64 _size;                   // total size of file
		uint64 _filled;                 // in fact required size of file

		uint64 _smers;                  // number of s-mers
		uint64 _kmers;                  // number of k-mers
		uint64 _ukmers;                 // number of u-mers

		uint64 _numberOfRuns;           // number of runs

	public:

		TempFile();

		void loadStats(std::string path, FILE *file);

		~TempFile();

		bool openW(const std::string filename);

		bool openR();

		bool write(SuperBundle *superBundle);

		bool write(char *data, const uint64 &size, const uint64 &smers, const uint64 &kmers, const uint64 &filled);

		bool read(SuperBundle *superBundle);

		uint64 approximateUniqueKmers(const double ratio) const;

		const uint64 &getKMersNumber() const;

		const uint64 &getSMersNumber() const;

		const uint64 &getUKMersNumber() const;

		const uint64 &getSize() const;

		inline const void calcNumberOfRuns(const double ratio, const uint64 maxUkmers);
		inline const void initNumberOfRuns() { _numberOfRuns = 1; }

		const uint64 &getNumberOfRuns() const;

		void incUKMersNumber(const uint64 &v);

		bool isEmpty();

		bool remove();

		void reset();

		void close();

		void fprintStat(FILE *file) const;
	};

	inline uint64 TempFile::approximateUniqueKmers(const double ratio) const {
		return _kmers * ratio;
	}

	inline const uint64 &TempFile::getKMersNumber() const {
		return _kmers;
	}

	inline const uint64 &TempFile::getSMersNumber() const {
		return _smers;
	}

	inline const uint64 &TempFile::getUKMersNumber() const {
		return _ukmers;
	}

	inline const uint64 &TempFile::getSize() const {
		return _size;
	}

	inline const uint64 &TempFile::getNumberOfRuns() const {
		return _numberOfRuns;
	}

	inline const void TempFile::calcNumberOfRuns(const double ratio, const uint64 maxUkmers) {
		uint64 x = approximateUniqueKmers(ratio);
		//_numberOfRuns = 1;
		_numberOfRuns = (x / (maxUkmers + 1)) + 1;
		//printf("calcNumberOfRuns: ratio=%f  maxUkMers= %lu  approximateUniqueKmers=%lu  _numberOfRuns=%lu\n", ratio, maxUkmers, x, _numberOfRuns);
	}

	inline void TempFile::incUKMersNumber(const uint64 &v) {
		__sync_add_and_fetch(&_ukmers, v);
	}

}

#endif /* BINFILE_H_ */
