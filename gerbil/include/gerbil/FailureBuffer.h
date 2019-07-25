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

#ifndef FAILUREBUFFER_H_
#define FAILUREBUFFER_H_

#include "Bundle.h"

namespace gerbil {
namespace cpu {

/**
 * Failure Buffer for cpu kmer bundles.
 * unsafe (single thread only)
 */
template<unsigned K>
class FailureBuffer {
	const size_t _maxBufferSize;	// max size of buffer in number of KMerBundle
	KMerBundle<K>** _buffer;		// buffer
	KMerBundle<K>** _top;			// top, next free item in buffer
	KMerBundle<K>** _end;			// end of buffer
	KMerBundle<K>* _currentBundle;	// current bundle
	FILE* _file;					// outsourcing file
	uint64 _amount;					// amount
	std::string _filePath;			// file path for outsourcing file

	typedef enum {fs_close, fs_rom, fs_wom} TFileState;
	TFileState _fileState;
	uint64 _storedNumber;

	void storeCurrentBundleToDisk();
	bool loadCurrentBundleFromDisk();
public:
	FailureBuffer(const size_t &maxBufferSize, std::string pPath, const uint &pId, const uint &pNr);
	~FailureBuffer();

	bool isEmpty();

	void clear();
	inline uint64 getAmount() const;

	// store
	void addKMer(const KMer<K> &kMer);

	// read
	bool getNextKMerBundle(KMerBundle<K>* &kMerBundle);
};

template<unsigned K>
void FailureBuffer<K>::storeCurrentBundleToDisk() {
	if(_fileState == fs_close) {
		_file = fopen (_filePath.c_str() , "wb+" );
		setbuf(_file, NULL);
		_fileState = fs_wom;
	}
	assert(_fileState == fs_wom);
	_currentBundle->store(_file);
	_currentBundle->clear();
	++_storedNumber;
}

template<unsigned K>
bool FailureBuffer<K>::loadCurrentBundleFromDisk() {
	if(!_storedNumber)
		return false;
	if(_fileState == fs_wom) {
		fclose(_file);
		_file = fopen (_filePath.c_str() , "rb" );
		setbuf(_file, NULL);
		_fileState = fs_rom;
	}
	assert(_fileState == fs_rom);
	// load stored bundle
	_currentBundle->load(_file);
	--_storedNumber;
	return true;
}

template<unsigned K>
FailureBuffer<K>::FailureBuffer(const size_t &maxBufferSize, std::string pPath, const uint &pId, const uint &pNr)
	:_maxBufferSize(maxBufferSize){
	_buffer = new KMerBundle<K>*[_maxBufferSize];
	_top = _buffer;
	_end = _buffer + _maxBufferSize;
	for(KMerBundle<K>** p(_buffer); p < _end; ++p)
		*p = new KMerBundle<K>();
	_currentBundle = new KMerBundle<K>();
	_file = 0;
	_filePath = pPath + "fails" + std::to_string(pId) + "_" + std::to_string(pNr);
	_fileState = fs_close;
	clear();
}

template<unsigned K>
FailureBuffer<K>::~FailureBuffer() {
	for(KMerBundle<K>** p(_buffer); p < _end; ++p)
			delete *p;
	delete[] _buffer;
}

template<unsigned K>
inline bool FailureBuffer<K>::isEmpty() {
	return (_top == _buffer) && _currentBundle->isEmpty();
}

template<unsigned K>
inline void FailureBuffer<K>::clear() {
	_amount = 0;
	_storedNumber = 0;
	///////////////
	while(_top > _buffer)
		(*(--_top))->clear();
	_currentBundle->clear();
	///////////////
	if(_fileState != fs_close) {
		fclose(_file);
		std::remove(_filePath.c_str());
		_fileState = fs_close;
	}
}

template<unsigned K>
inline uint64 FailureBuffer<K>::getAmount() const {
	return _amount;
}

// save kMer
template<unsigned K>
void FailureBuffer<K>::addKMer(const KMer<K> &kMer) {
	if(!_currentBundle->add(kMer)) {
		if(_top != _end) {	// save in buffer
			std::swap(*_top, _currentBundle);
			++_top;
		}
		else	// save on disk
			storeCurrentBundleToDisk();
		_currentBundle->add(kMer);
	}
	++_amount;
}

template<unsigned K>
bool FailureBuffer<K>::getNextKMerBundle(KMerBundle<K>* &kMerBundle) {
	if(_top != _buffer) { //from buffer
		kMerBundle = *(--_top);
		return true;
	}
	if(_currentBundle->isEmpty()) // load from disk
		if(!loadCurrentBundleFromDisk())
			return false;
	kMerBundle = _currentBundle;
	//printf("#");
	return true;

}

}


/**
 * Failure Buffer for gpu kmer bundles.
 */

namespace gpu {

// buffer for failure kMers
// unsafe (single thread only)
template<unsigned K>
class FailureBuffer {
	const size_t _maxBufferSize;	// max size of buffer in number of KMerBundle
	KMerBundle<K>** _buffer;		// buffer
	KMerBundle<K>** _top;			// top, next free item in buffer
	KMerBundle<K>** _end;			// end of buffer
	KMerBundle<K>* _currentBundle;	// current bundle
	FILE* _file;
	uint64 _amount;
	std::string _filePath;

	typedef enum {fs_close, fs_rom, fs_wom} TFileState;
	TFileState _fileState;
	uint64 _storedNumber;

	void storeCurrentBundleToDisk();
	bool loadCurrentBundleFromDisk();
public:
	FailureBuffer(const size_t &maxBufferSize, std::string pPath, const uint &pId, const uint &pNr);
	~FailureBuffer();

	bool isEmpty();

	void clear();
	inline uint64 getAmount() const;

	// store
	void addKMer(const KMer<K> &kMer);

	// read
	bool getNextKMerBundle(KMerBundle<K>* &kMerBundle);
};

template<unsigned K>
void FailureBuffer<K>::storeCurrentBundleToDisk() {
	//printf("-------store to disk\n");
	if(_fileState == fs_close) {
		_file = fopen (_filePath.c_str() , "wb+" );
		setbuf(_file, NULL);
		_fileState = fs_wom;
	}
	assert(_fileState == fs_wom);
	_currentBundle->store(_file);
	_currentBundle->clear();
	++_storedNumber;
}

template<unsigned K>
bool FailureBuffer<K>::loadCurrentBundleFromDisk() {
	if(!_storedNumber)
		return false;
	//printf("-------load from disk\n");
	if(_fileState == fs_wom) {
		fclose(_file);
		_file = fopen (_filePath.c_str() , "rb" );
		setbuf(_file, NULL);
		_fileState = fs_rom;
	}
	assert(_fileState == fs_rom);
	// load stored bundle
	_currentBundle->load(_file);
	--_storedNumber;
	return true;
}

template<unsigned K>
FailureBuffer<K>::FailureBuffer(const size_t &maxBufferSize, std::string pPath, const uint &pId, const uint &pNr)
	:_maxBufferSize(maxBufferSize){
	_top = _buffer = new KMerBundle<K>*[_maxBufferSize];
	_end = _buffer + _maxBufferSize;
	for(KMerBundle<K>** p(_buffer); p != _end; ++p)
		*p = new KMerBundle<K>();
	_currentBundle = new KMerBundle<K>;
	_file = 0;
	_filePath = pPath + "fails" + std::to_string(pId) + "_" + std::to_string(pNr);
	_fileState = fs_close;
	clear();
}

template<unsigned K>
FailureBuffer<K>::~FailureBuffer() {
	for(KMerBundle<K>** p(_buffer); p < _end; ++p)
		delete *p;
	delete[] _buffer;
	//delete _currentBundle;
}

template<unsigned K>
inline bool FailureBuffer<K>::isEmpty() {
	return (_top == _buffer) && _currentBundle->isEmpty();
}

template<unsigned K>
inline void FailureBuffer<K>::clear() {
	_amount = 0;
	_storedNumber = 0;
	///////////////
	while(_top > _buffer)
		(*(--_top))->clear();
	_currentBundle->clear();
	///////////////
	if(_fileState != fs_close) {
		fclose(_file);
		std::remove(_filePath.c_str());
		_fileState = fs_close;
	}
}

template<unsigned K>
inline uint64 FailureBuffer<K>::getAmount() const {
	return _amount;
}

// save kMer
template<unsigned K>
void FailureBuffer<K>::addKMer(const KMer<K> &kMer) {
	if(!_currentBundle->add(kMer)) {
		if(_top != _end) {	// save in buffer
			std::swap(*_top, _currentBundle);
			++_top;
		}
		else	// save on disk
			storeCurrentBundleToDisk();
		_currentBundle->add(kMer);
	}
	++_amount;
}

template<unsigned K>
bool FailureBuffer<K>::getNextKMerBundle(KMerBundle<K>* &kMerBundle) {
	if(_top != _buffer) { //from buffer
		kMerBundle = *(--_top);
		return true;
	}
	if(_currentBundle->isEmpty()) // load from disk
		if(!loadCurrentBundleFromDisk())
			return false;
	kMerBundle = _currentBundle;
	//printf("#");
	return true;

}

}


}

#endif /* FAILUREBUFFER_H_ */
