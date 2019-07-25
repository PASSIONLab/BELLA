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

#ifndef SUPERWRITER_HPP_
#define SUPERWRITER_HPP_

#include "SyncQueue.h"
#include "TempFile.h"

namespace gerbil {

/*
 * writes SuperBundles to disk
 */
class SuperWriter {
	SyncSwapQueueMPSC<SuperBundle>** _superBundleQueues;

	std::thread** _processThreads;	// number of worker threads

	TempFile* _tempFiles;			// number of temporary files
	uint_tfn _tempFilesNumber;		// number of temporary files

	uint64 _maxBufferSize;			// maximum buffer size

	uint64 _superBundlesNumber;		// number of SuperBundles
	uint64 _tempFilesFilledSize;	// real use of space
	uint64 _sMersNumber;			// number of s-mers
	uint64 _kMersNumber;			// number of k-mers

public:
	TempFile* getTempFiles();	// returns the temporary files

	/*
	 * constructor
	 */
	SuperWriter(std::string pTempFolder,
			SyncSwapQueueMPSC<SuperBundle>** superBundleQueues,
			const uint_tfn &tempFilesNumber, const uint64 &maxBufferSize);

    /*
     * starts the entire working process
     */
	void process();

  	/*
   	 * joins all threads
   	 */
	void join();

	/*
     * prints some statistical outputs
     */
	void print();

	/*
	 * destructor
	 */
	~SuperWriter();
};


////////////////////////////////////////////////////////////////////////
// beginning of the test region
////////////////////////////////////////////////////////////////////////

//#define SWB_HEAP

class SuperBundleStackItem {
public:
	SuperBundle* superBundle;
	SuperBundleStackItem* next;
	SuperBundleStackItem(SuperBundle* sb) :
			superBundle(sb), next(NULL) {
	}
	~SuperBundleStackItem() {
		delete superBundle;
	}
};

class SuperBundleStack {
	SuperBundleStackItem* _top;
	uint64 _size;
public:
#ifdef SWB_HEAP
	SuperBundleStack* parent;
	SuperBundleStack* leftChild;
	SuperBundleStack* rightChild;
#endif
	SuperBundleStack() :
			_top(NULL), _size(0)
#ifdef SWB_HEAP
	,parent(NULL)
	,leftChild(NULL)
	,rightChild(NULL)
#endif
	{
	}

	inline const uint64& size() const {
		return _size;
	}

	inline void push(SuperBundleStackItem* &superBundleStackItem) {
		superBundleStackItem->next = _top;
		_top = superBundleStackItem;
		++_size;
	}

	inline void push(SuperBundle* superBundle) {
		SuperBundleStackItem* superBundleStackItem = new SuperBundleStackItem(
				superBundle);
		push(superBundleStackItem);
	}

	inline bool pop(SuperBundleStackItem* &superBundleStackItem) {
		if (_top) {
			superBundleStackItem = _top;
			_top = _top->next;
			superBundleStackItem->next = NULL;
			--_size;
			return true;
		}
		return false;
	}

	inline void swap(SuperBundleStack* &other) {
		uint64 s(other->_size);
		other->_size = _size;
		_size = s;
		SuperBundleStackItem* sbsi(other->_top);
		other->_top = _top;
		_top = sbsi;
	}
};

class SuperBundleConcatenator {
	SuperBundleStack* superBundleStacks;		// by id of Tempfiles
	SuperBundleStack emptySuperBundleStack;		// stack for empty SuperBundles
	const uint_tfn _tempFilesNumber;			// number of TempFiles
	uint64 _maxEmptySuperBundlesNumber;		// max number of empty SuperBundles
	uint64 _emptySuperBundlesNumber;

#ifdef SWB_HEAP
	SuperBundleStack* _root;

	void setupHeap(SuperBundleStack* parent, const uint64 &pId) {
		uint64 childId = (pId << 1) | 1;
		if(childId < _tempFilesNumber) {
			parent->leftChild = superBundleStacks + childId;
			parent->leftChild->parent = parent;
			setupHeap(parent->leftChild, childId);
		}
		if(++childId < _tempFilesNumber) {
			parent->rightChild = superBundleStacks + childId;
			parent->rightChild->parent = parent;
			setupHeap(parent->rightChild, childId);
		}
	}
#else
	uint64 _lastStackNr;
#endif
public:
	SuperBundleConcatenator(const uint_tfn &tempFilesNumber,
			uint64 maxEmptySuperBundlesNumber) :
			_tempFilesNumber(tempFilesNumber), _maxEmptySuperBundlesNumber(
					maxEmptySuperBundlesNumber), _emptySuperBundlesNumber(0) {
		superBundleStacks = new SuperBundleStack[_tempFilesNumber];

#ifdef SWB_HEAP
		_root = superBundleStacks;
		setupHeap(_root, 0);
#else
		_lastStackNr = 0;
#endif

		uint64 t = std::min(_tempFilesNumber << 1, _maxEmptySuperBundlesNumber >> 4);
		do {
			emptySuperBundleStack.push(new SuperBundle);
			++_emptySuperBundlesNumber;
		} while (_emptySuperBundlesNumber < t);
	}
	~SuperBundleConcatenator() {
		std::printf("||||||||||||||||||||||%7lu / %7lu\n",
				_emptySuperBundlesNumber, _maxEmptySuperBundlesNumber);

		SuperBundleStackItem* sbsi;
		while (emptySuperBundleStack.size()) {
			emptySuperBundleStack.pop(sbsi);
			delete sbsi;
		}

		delete[] superBundleStacks;
	}

	inline bool notFull() {
		if (emptySuperBundleStack.size())
			return true;
		else if (_emptySuperBundlesNumber >= _maxEmptySuperBundlesNumber)
			return false;
		uint64 t = std::min(_emptySuperBundlesNumber + _tempFilesNumber,
				_maxEmptySuperBundlesNumber);
		do {
			emptySuperBundleStack.push(new SuperBundle);
			++_emptySuperBundlesNumber;
		} while (_emptySuperBundlesNumber < t);
		return true;
	}

	// with respect to _root
	inline bool notEmpty() {
#ifdef SWB_HEAP
		return _root->size();
#else
		return _emptySuperBundlesNumber != emptySuperBundleStack.size();
#endif
	}

	// with respect to _root
	inline bool isEmpty() {
		return !notEmpty();
	}

	// Cond.: notFull()
	bool swapPush(SuperBundle* &superBundle) {
		SuperBundleStackItem* sbsi;
		if (!emptySuperBundleStack.pop(sbsi))
			return false;
		uint64 tid = superBundle->tempFileId;

		// swap SuperBundle (full <--> empty)
		SuperBundle* sb = sbsi->superBundle;
		sbsi->superBundle = superBundle;
		superBundle = sb;

		SuperBundleStack* sbs = superBundleStacks + tid;

		// add SuperBundle to Stack
		sbs->push(sbsi);

#ifdef SWB_HEAP
		SuperBundleStack* psbs = sbs->parent;

		// update heap
		while(psbs && sbs->size() > psbs->size()) {
			sbs->swap(psbs);
			sbs = psbs;
			psbs = psbs->parent;
		}
#endif

		return true;
	}

	// Cond.: notEmpty()
	bool swapPop(SuperBundle* &superBundle) {
		if (isEmpty())
			return false;

		SuperBundleStackItem* sbsi;

		// pop item from stack
#ifdef SWB_HEAP
		_root->pop(sbsi);
#else
		SuperBundleStack* sbs;
		while (!(sbs = superBundleStacks + (_lastStackNr++ % _tempFilesNumber))->size())
			;
		sbs->pop(sbsi);
#endif

		// swap SuperBundle (empty <--> full)
		SuperBundle* sb = sbsi->superBundle;
		sbsi->superBundle = superBundle;
		superBundle = sb;

		// push empty superBundle to stack
		emptySuperBundleStack.push(sbsi);

		return true;
	}

#ifdef SWB_HEAP
	// update key for root (after each swapPop series)
	void updateHeap() {
		SuperBundleStack* sbs = _root;
		SuperBundleStack* lcsbs = sbs->leftChild;
		SuperBundleStack* rcsbs = sbs->rightChild;

		// update
		do {
			// left child greater
			if(lcsbs && lcsbs->size() > sbs->size()) {
				sbs->swap(lcsbs);
				sbs = lcsbs;
				rcsbs = lcsbs->rightChild;
				lcsbs = lcsbs->leftChild;
			}
			// right child greater
			else if(rcsbs && rcsbs->size() > sbs->size()) {
				sbs->swap(rcsbs);
				sbs = rcsbs;
				lcsbs = rcsbs->leftChild;
				rcsbs = rcsbs->rightChild;
			}
			else
			break;
		}while(true);
	}
#endif

	void print() {
#ifdef SWB_HEAP
		printf("HEAP:\n");
		_print(_root, 0);
#endif
	}
private:
#ifdef SWB_HEAP
	void _print(SuperBundleStack* node, uint64 depth) {
		if(!node)
		return;
		for(uint64 i = 0; i < depth; ++i)
		printf("\t");
		printf("[%5lu]\n", node->size());
		++depth;
		_print(node->leftChild, depth);
		_print(node->rightChild, depth);

	}
#endif
};

}

#endif /* SUPERWRITER_HPP_ */
