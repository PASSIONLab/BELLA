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

#ifndef FASTREADER_H_
#define FASTREADER_H_

#include "global.h"
#include "SyncQueue.h"
#include "Bundle.h"
#include "FastFile.h"

#include <fstream>
#include <zlib.h>
#include <bzlib.h>

namespace gerbil {

/*
 * reads files from disk and stored the data as byte blocks
 * decompresses the files if necessary
 */
	class FastReader {

	private:
		uint64 _totalReadBytes;                            // total number of bytes which have been read
		bfs::path _path;                                   // file path
		TFileType _fileType;                               // file type

		FastFile **_fastFiles;                             // list of FastFiles
		uint_fast32_t _fastFilesNumber;                    // number of FastFiles

		SyncSwapQueueSPSC<FastBundle> **_syncSwapQueues;   // SyncSwapQueues for FastBundles

		std::thread **_processThreads;                     // worker threads

		std::atomic<uint64> _totalBlocksRead;              // total number of blocks which have been read

		uint8 _threadsNumber;                              // number of threads

		std::atomic<uint64> _fastFileNr;                   // number of recently processed file

		/*
		 * reads a single file and pushes the FastBundles in the SyncSwapQueue
		 */
		void readFile(
				const uint tId,
				const FastFile &fastFile,
                SyncSwapQueueSPSC<FastBundle> &syncSwapQueue
#ifdef DEB_MESS_FASTREADER
				, StopWatch* sw
#endif
		);

		/*
		 * starts the working process of a single thread
		 */
		void processThread(const size_t tId,
		                   SyncSwapQueueSPSC<FastBundle> &syncSwapQueue);

	public:
		/*
		 * constructor
		 */
		FastReader(const uint32 &frBlocksNumber, std::string pPath,
		           uint8 &_readerParserThreadsNumber);

		SyncSwapQueueSPSC<FastBundle> **getSyncSwapQueues(); // returns the list of SyncSwapQueues
		TFileType getFileType() const;                       // returns the file type

		/*
		 * starts the entire working process
		 */
		void process();

		/*
		 * prints some statistical outputs
		 */
		void print();

		/*
		 * joins all threads
		 */
		void join();

		/*
		 * destructor
		 */
		~FastReader();
	};

}

#endif /* FASTREADER_H_ */
