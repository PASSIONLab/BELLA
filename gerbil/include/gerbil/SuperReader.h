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

#ifndef SUPERREADER_H_
#define SUPERREADER_H_

#include "Bundle.h"
#include "SyncQueue.h"
#include "TempFile.h"
#include "KmerDistributer.h"

namespace gerbil {

/*
 * reads SuperBundles from disk
 */
class SuperReader {
	SyncSwapQueueSPMC<SuperBundle> _syncSwapQueue;	// SuperBundleQueue

	TempFile* _tempFiles;							// temporary files
	uint_tfn _tempFilesNumber;						// number of temporary files
	uint_tfn* _tempFilesOrder;						// processing order of temporary files

	std::thread* _processThread;					// process thread

	KmerDistributer* _distributor;
	
	/*
	 * starts the working process of a single thread
	 */
	void processThread();
public:

	uint_tfn* getTempFilesOrder();							// returns the order of temporary files
	SyncSwapQueueSPMC<SuperBundle>* getSuperBundleQueue();	// returns the SuperBundleQueue

	/*
	 * constructor
	 */
	SuperReader(
			const uint32 &superBundlesNumber,
			TempFile* tempFiles,
			const uint_tfn& tempFilesNumber,
			KmerDistributer* distributor
	);

    /*
     * starts the entire working process
     */
	void process();

  	/*
   	 * joins all threads
   	 */
	void join();

	/*
	 * constructor
	 */
	~SuperReader();
};

}


#endif /* SUPERREADER_H_ */
