/*
 * SuperWriter.cpp
 *
 *  Created on: 20.05.2015
 *      Author: marius
 */

#include "../../include/gerbil/SuperWriter.h"

gerbil::SuperWriter::SuperWriter(std::string pTempFolder,
		SyncSwapQueueMPSC<SuperBundle>** superBundleQueues,
		const uint_tfn &tempFilesNumber, const uint64 &maxBufferSize) :
		_superBundleQueues(superBundleQueues), _tempFilesNumber(
				tempFilesNumber), _maxBufferSize(maxBufferSize), _processThreads(
		NULL), _superBundlesNumber(0), _tempFilesFilledSize(0), _sMersNumber(0), _kMersNumber(
				0) {
	// open all temp files
	std::string tempPathName;
	_tempFiles = new TempFile[_tempFilesNumber];
	for (uint_tfn tempFileId = 0; tempFileId < _tempFilesNumber; ++tempFileId) {
		tempPathName = pTempFolder + "temp" + std::to_string(tempFileId)
				+ ".bin";
		if (!_tempFiles[tempFileId].openW(tempPathName)) {
			std::cerr << "unable to create temp File '" << tempPathName << "'"
					<< std::endl;
			exit(1);
		}
	}

	_processThreads = new std::thread*[SB_WRITER_THREADS_NUMBER];
}

gerbil::SuperWriter::~SuperWriter() {
	//delete[] _tempFiles;
	// _tempFiles is been deleted in SuperReader.cpp
}

#define storeSB(sb) {												\
	_tempFilesFilledSize += sb->getSize();							\
	_kMersNumber += sb->kMerNumber;									\
	_sMersNumber += sb->sMerNumber;									\
	++_superBundlesNumber;											\
	if(!_tempFiles[sb->tempFileId].write(sb)) {						\
		std::cerr << "write super-mers to temp-file failed" << std::endl;		\
		exit(74);													\
	}																\
	sb->clear();													\
}

#define storeSB2(sb) {												\
	tempFilesFilledSize += sb->getSize();							\
	kMersNumber += sb->kMerNumber;									\
	sMersNumber += sb->sMerNumber;									\
	++superBundlesNumber;											\
	if(!_tempFiles[sb->tempFileId].write(sb)) {						\
		std::cerr << "write super-mers to temp-file failed" << std::endl;		\
		exit(74);													\
	}																\
	sb->clear();													\
}

void gerbil::SuperWriter::process() {
	for (uint i = 0; i < SB_WRITER_THREADS_NUMBER; i++)
		_processThreads[i] =
				new std::thread([this](uint64 tId) {
					IF_MESS_SUPERWRITER(
							StopWatch sw;
							sw.start();
					)

					uint64 n = 0; uint64 s = 0;

// best first
#if false
						{
							printf("------------->_maxBufferSize: %7lu\n", _maxBufferSize);
							SuperBundleConcatenator superBundleConcatenator(_tempFilesNumber, _maxBufferSize);
							SuperBundle* sb = new SuperBundle;
							while(true) {
								// fill SuperBundleConcatenator
								while(superBundleConcatenator.notFull() && _superBundleQueue->swapPop_nl(sb))
								superBundleConcatenator.swapPush(sb);
								// wait for at least one item
								if(superBundleConcatenator.isEmpty()) {
									IF_MESS_SUPERWRITER(sw.hold();)
									if(_superBundleQueue->swapPop(sb)) {
										IF_MESS_SUPERWRITER(sw.proceed();)
										superBundleConcatenator.swapPush(sb);
									}
									else {
										IF_MESS_SUPERWRITER(sw.proceed();)
										break;
									}
								}
								// store 1 stack
								while(superBundleConcatenator.notEmpty()) {
									superBundleConcatenator.swapPop(sb);
									storeSB(sb);
									++s;
								}
								++n;
#ifdef SWB_HEAP
						superBundleConcatenator.updateHeap();
#endif
					}
					delete sb;
				}
#endif

// with greedy
#if false
						SuperBundle* esb;
						SuperBundle* emptySuperBundles[_maxBufferSize];
						SuperBundle** emptySuperBundles_last = emptySuperBundles;

						for(uint64 i(0); i < _maxBufferSize; ++i)
						*(emptySuperBundles_last++) = new SuperBundle();

						esb = *(--emptySuperBundles_last);

						queue<SuperBundle*> superBundleBuffer[_tempFilesNumber];

						size_t p = 0;

						uint64 c(0);

						while(true/*!_superBundleQueue->isFinalized()*/) {
							while(emptySuperBundles_last > emptySuperBundles && _superBundleQueue->swapPop_nl(esb)) {
								superBundleBuffer[esb->tempFileId].push(esb);
								esb = *(--emptySuperBundles_last);
								++c;
								//printf("push\n");
							}
							if(!c && emptySuperBundles_last > emptySuperBundles) {
								IF_MESS_SUPERWRITER(sw.hold();)
								if(_superBundleQueue->swapPop(esb)) {
									IF_MESS_SUPERWRITER(sw.proceed();)
									//printf("wait...\n");
									superBundleBuffer[esb->tempFileId].push(esb);
									esb = *(--emptySuperBundles_last);
									++c;
									//printf("push\n");
									//printf("wait res=%lu\n", c);
								}
								else {
									IF_MESS_SUPERWRITER(sw.proceed());
									break;
								}
							}
							if(c) {
								queue<SuperBundle*>* sbbq;
								while((sbbq = superBundleBuffer + (p++ % _tempFilesNumber))->empty());
								//uint64 x = 0;
								while(!sbbq->empty()) {
									//++x;
									SuperBundle* sb = sbbq->front();
									sbbq->pop();
									storeSB(sb);
									*(emptySuperBundles_last++) = sb;
									--c;
									//printf("pop\n");
									++s;
								}
								++n;
								//if(x) printf("%4lu --> %4lu / %4lu   : %7lu", (p-1)% _tempFilesNumber, x, c, _superBundlesNumber);
								//putchar('\n');
							}
						}
						//printf("finalized!\n");
						/*for(size_t i(0); i < _tempFilesNumber; ++i) {
						 queue<SuperBundle*>* sbbq = superBundleBuffer + i;
						 while(!sbbq->empty()) {
						 SuperBundle* sb = sbbq->front();
						 sbbq->pop();
						 storeSB(sb);
						 --c;
						 ++s;
						 delete sb;
						 }
						 ++n;
						 }*/

						delete esb;
						while(emptySuperBundles_last > emptySuperBundles) {
							esb = *(--emptySuperBundles_last);
							delete esb;
						}
#endif

// default
#if false
						SuperBundle* sb = new SuperBundle;
						// swappop all superbundles from queue
						IF_MESS_SUPERWRITER(sw.hold();)
						while(_superBundleQueue->swapPop(sb)) {
							IF_MESS_SUPERWRITER(sw.proceed();)
							storeSB(sb);
							++n;
							++s;
							IF_MESS_SUPERWRITER(sw.hold();)
						}
						IF_MESS_SUPERWRITER(sw.proceed();)
						delete sb;
#endif

						// multifile
#if true

						uint64 tempFilesFilledSize = 0;
						uint64 kMersNumber = 0;
						uint64 sMersNumber = 0;
						uint64 superBundlesNumber = 0;
						SuperBundle* sb = new SuperBundle;
						// swappop all superbundles from queue
						IF_MESS_SUPERWRITER(sw.hold();)
						SyncSwapQueueMPSC<SuperBundle>& superBundleQueue = *_superBundleQueues[tId];

						// greedy buffering
#if false
						SuperBundle* esb;
						SuperBundle* emptySuperBundles[_maxBufferSize];
						SuperBundle** emptySuperBundles_last = emptySuperBundles;

						// max buffer size for a single write operation
						const uint64 maxPartBufferSize = 64 * 1024 * 1024;// 64MB
						char* buffer;
						buffer = new char[maxPartBufferSize];

						// initialize empty buffer
						for(uint64 i(0); i < _maxBufferSize; ++i)
						*(emptySuperBundles_last++) = new SuperBundle();

						esb = *(--emptySuperBundles_last);

						// buffer queues
						std::queue<SuperBundle*> superBundleBuffer[_tempFilesNumber];

						size_t p = 0;

						uint64 c(0);
						while(true) {
							// transfer superbundles from main queue to buffer queues
							while(emptySuperBundles_last > emptySuperBundles && superBundleQueue.swapPop_nl(esb)) {
								superBundleBuffer[esb->tempFileId].push(esb);
								esb = *(--emptySuperBundles_last);
								++c;
							}
							// wait for at least 1 superbundle (or break loop)
							if(!c && emptySuperBundles_last > emptySuperBundles) {
								IF_MESS_SUPERWRITER(sw.hold();)
								if(superBundleQueue.swapPop(esb)) {
									IF_MESS_SUPERWRITER(sw.proceed();)
									//printf("wait...\n");
									superBundleBuffer[esb->tempFileId].push(esb);
									esb = *(--emptySuperBundles_last);
									++c;
								}
								else {
									IF_MESS_SUPERWRITER(sw.proceed());
									break;
								}
							}
							if(c) {
								std::queue<SuperBundle*>* sbbq;
								while((sbbq = superBundleBuffer + (p++ % _tempFilesNumber))->empty());
								uint64 curBufferSize = 0;
								uint64 smers = 0;
								uint64 kmers = 0;
								uint64 fill = 0;
								uint64 tempId = sbbq->front()->tempFileId;
								while(!sbbq->empty() && curBufferSize + SUPER_BUNDLE_DATA_SIZE_B <= maxPartBufferSize) {
									//++x;
									SuperBundle* sb = sbbq->front();
									sbbq->pop();

									tempFilesFilledSize += sb->getSize();
									kmers += sb->kMerNumber;
									smers += sb->sMerNumber;
									fill += sb->getSize();
									++superBundlesNumber;
									memcpy(buffer + curBufferSize, sb->data,SUPER_BUNDLE_DATA_SIZE_B);
									curBufferSize += SUPER_BUNDLE_DATA_SIZE_B;
									sb->clear();

									*(emptySuperBundles_last++) = sb;
									--c;
									++s;
								}
								kMersNumber += kmers;
								sMersNumber += smers;
								if(!_tempFiles[tempId].write(buffer, curBufferSize, smers, kmers, fill)) {
									std::cerr << "write buffered super-mers to temp-file failed" << std::endl;
									exit(75);
								}

								++n;
							}

						}

						delete esb;
						while(emptySuperBundles_last > emptySuperBundles) {
							esb = *(--emptySuperBundles_last);
							delete esb;
						}

						delete[] buffer;

#endif

						// default
#if true
						while(superBundleQueue.swapPop(sb)) {
							IF_MESS_SUPERWRITER(sw.proceed();)
							storeSB2(sb);
							++n;
							++s;
							IF_MESS_SUPERWRITER(sw.hold();)
						}
#endif

						__sync_add_and_fetch(&_tempFilesFilledSize, tempFilesFilledSize);
						__sync_add_and_fetch(&_kMersNumber, kMersNumber);
						__sync_add_and_fetch(&_sMersNumber, sMersNumber);
						__sync_add_and_fetch(&_superBundlesNumber, superBundlesNumber);

						IF_MESS_SUPERWRITER(sw.proceed();)
						delete sb;
#endif

						//printf(">>>>>%5.2f\n", (double)s / n);

						IF_DEB(printf("close files\n"));
						// close all temp files
						for(uint_tfn tempFileId(0); tempFileId < _tempFilesNumber; ++tempFileId)
						if(tempFileId % SB_WRITER_THREADS_NUMBER == tId)
						_tempFiles[tempFileId].close();
						IF_MESS_SUPERWRITER(
								sw.stop();
								printf("time SuperWriter: %.3f s\n", sw.get_s());
						)
					}, i);
}

gerbil::TempFile* gerbil::SuperWriter::getTempFiles() {
	return _tempFiles;
}

void gerbil::SuperWriter::join() {
	for (uint i = 0; i < SB_WRITER_THREADS_NUMBER; ++i) {
		_processThreads[i]->join();
		delete _processThreads[i];
	}IF_DEB(printf("superWriter is rdy...\n"));
}

void gerbil::SuperWriter::print() {
	printf("number of superbundles : %12lu\n", _superBundlesNumber);
	printf("size of temp file      : % 12.3f MB (%5.3f MB, %3.1f%% belegt)\n",
			(double) _superBundlesNumber * (SUPER_BUNDLE_DATA_SIZE_B / 1024)
					/ 1024, (double) _tempFilesFilledSize / 1024 / 1024,
			((double) _tempFilesFilledSize)
					/ ((double) SUPER_BUNDLE_DATA_SIZE_B * _superBundlesNumber)
					* 100);
	printf("number of s-mers       : %12lu\n", _sMersNumber);
	printf("number of k-mers       : %12lu\n", _kMersNumber);
}

