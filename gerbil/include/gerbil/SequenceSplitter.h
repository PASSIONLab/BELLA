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

#ifndef SUPERSPLITTER_H_
#define SUPERSPLITTER_H_

#include "SyncQueue.h"
#include "Bundle.h"

namespace gerbil {

/*
 * splits reads of ReadBundles to s-mers and stores them into several SuperBundleQueues
 */
class SequenceSplitter{
	SyncSwapQueueMPMC<ReadBundle>* _readBundleSyncQueue;   // ReadBundleQueue
	SyncSwapQueueMPSC<SuperBundle>** _superBundleQueues;   // SuperBundleQueues

    byte _splitterThreadCount;         // number of splitter threads
    std::thread** _splitterThreads;    // list of splitter threads

	uint32_t _k;                       // size of k-mer
	byte _m;                           // size of minimizer
	uint_tfn _tempFilesNumber;         // number of temporary files

	uint32* _mVal;                     // minimizer values for the ranking
	uint_tfn* _mToBin;                 // assignment of minimizer to temporary files
	const bool _norm;                  // checked if the normalization is enabled

	SuperBundle*** _leftSuperBundles;  // incomplete SuperBundles

	std::atomic<uint64> _baseNumbers;  // total number of bases

	uint32_t invMMer(const uint32_t &mmer);    // inverts a minimizer
	bool isAllowed(uint32 mmer);               // checks whether a minimizer is allowed (special, tested strategies)
	void detMMerHisto();                       // calculation of a histogram (special, tested strategies)

	void processThread(const uint &id);        // starts a single thread

public:
    SyncSwapQueueMPSC<SuperBundle>** getSuperBundleQueues();    // returns the SuperBundleQueues

    /*
     * constructor
     */
	SequenceSplitter(
			const uint32 &superBundlesNumber,						// size of superBundleQueue
			SyncSwapQueueMPMC<ReadBundle>* readBundleSyncQueue,		// Queue with Bundles of reads
			const uint8 &splitterThreadsNumber,						// Number of SplitterThreads
			const uint32_t &k,										// size of k-mer
			const uint8 &m,											// Size of m-mer (minimizer)
			const uint_tfn &tempFilesNumber,						// Number of Bin-Files
			const bool &norm										// normalized k-mers
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
     * prints some statistical outputs
     */
    void print();

    /*
     * destructor
     */
    ~SequenceSplitter();
};


}


#endif /* SUPERSPLITTER_H_ */
