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

#ifndef TEMPFILESTATISTIC_H_
#define TEMPFILESTATISTIC_H_

#include <cmath>
#include "TempFile.h"

namespace gerbil {

class TempFileStatistic {
	const TempFile* const _tempFiles;
	const uint_tfn _tempFilesNumber;
	uint64 _minKMersNumber, _maxKMersNumber, _avgKMersNumber, _varKMersNumber, _sdKMersNumber, _sumKMersNumber;
	uint64 _maxSize,  _avgSize, _varSize, _sdSize, _sumSize;

	// tempFilesNumber > 2
	void calcStatistic() {
		_minKMersNumber = _maxKMersNumber = _avgKMersNumber = _sumKMersNumber = _tempFiles[0].getKMersNumber();
		_maxSize = _avgSize = _sumSize = _tempFiles[0].getSize();

		uint64 curKMersNumber, curSize;
		const TempFile* const tempFile_p_end = _tempFiles + _tempFilesNumber;
		for(const TempFile* tempFile_p(_tempFiles + 1); tempFile_p < tempFile_p_end; ++tempFile_p){
			curKMersNumber = tempFile_p->getKMersNumber();
			_sumKMersNumber += curKMersNumber;
			if(curKMersNumber > _maxKMersNumber)
				_maxKMersNumber = curKMersNumber;
			else if(curKMersNumber < _minKMersNumber)
				_minKMersNumber = curKMersNumber;

			curSize = tempFile_p->getSize();
			_sumSize += curSize;
			if(curSize > _maxSize)
				_maxSize = curSize;
		}
		_avgKMersNumber = _sumKMersNumber / _tempFilesNumber;
		_avgSize = _sumSize / _tempFilesNumber;

		_varKMersNumber = 0;
		_varSize = 0;
		int64 tfn_s1 = _tempFilesNumber - 1;
		int64 diffKMersNumber, diffSize;
		for(const TempFile* tempFile_p(_tempFiles); tempFile_p < tempFile_p_end; ++tempFile_p){
			diffKMersNumber = tempFile_p->getKMersNumber();
			diffKMersNumber -= _avgKMersNumber;
			_varKMersNumber += diffKMersNumber / tfn_s1 * diffKMersNumber;

			diffSize = tempFile_p->getSize();
			diffSize -= _avgSize;
			_varSize += diffSize / tfn_s1 * _varSize;
		}

		_sdKMersNumber = sqrt(_varKMersNumber);
		_sdSize = sqrt(_varSize);
	}
public:
	TempFileStatistic(TempFile* pTempFiles, uint_tfn pTempFilesNumber)
	: _tempFiles(pTempFiles), _tempFilesNumber(pTempFilesNumber){
		calcStatistic();
	}

	inline const uint64& getMinKMersNumber() const {
		return _minKMersNumber;
	}

	inline const uint64& getMaxKMersNumber() const {
		return _maxKMersNumber;
	}

	inline const uint64& getAvgKMersNumber() const {
		return _avgKMersNumber;
	}

	inline const uint64& getVarKMersNumber() const {
		return _varKMersNumber;
	}

	inline const uint64& getSdKMersNumber() const {
		return _sdKMersNumber;
	}

	inline const uint64& getSumKMersNumber() const {
		return _sumKMersNumber;
	}

	inline uint64 getAvg1SdKMersNumber() const {
		return _avgKMersNumber + _sdKMersNumber;
	}

	inline uint64 getAvg2SdKMersNumber() const {
		return _avgKMersNumber + 2 * _sdKMersNumber;
	}

	inline uint64 getAvg3SdKMersNumber() const {
		return _avgKMersNumber + 3 * _sdKMersNumber;
	}

	inline const uint64& getSumSize() const {
		return _sumSize;
	}

	inline uint64 getAvg2SdSize() const {
		return _avgSize + 2 * _sdSize;
	}
};

}

#endif /* TEMPFILESTATISTIC_H_ */
