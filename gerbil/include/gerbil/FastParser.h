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

#ifndef FASTPARSER_H_
#define FASTPARSER_H_

#include "SyncQueue.h"
#include "Bundle.h"
#include "global.h"

#include <zlib.h>
#include <bzlib.h>

namespace gerbil {

	class Decompressor {

	private:
		int ret;
		z_stream strm;
		FastBundle* _compressedFastBundle;
		bool _resetable;

	public:
		FastBundle fastBundle;

		void reset() {
			if(!_resetable) {
				_resetable = true;
				return;
			}
			//std::cout << "reset Decompressor" << std::endl;
			/* allocate inflate state */
			strm.zalloc = Z_NULL;
			strm.zfree = Z_NULL;
			strm.opaque = Z_NULL;
			strm.avail_in = 0;
			strm.next_in = Z_NULL;
			//ret = inflateReset(&strm);
			ret = inflateReset2(&strm, 16+MAX_WBITS);
			if (ret != Z_OK) {
				std::cerr << "decompress stream: unkown error" << std::endl;
				exit(1);
			}
			fastBundle.clear();
			fastBundle.compressType = fc_DECOMPRESSOR;
		}

		Decompressor() {
			//std::cout << "init Decompressor" << std::endl;
			/* allocate inflate state */
			strm.zalloc = Z_NULL;
			strm.zfree = Z_NULL;
			strm.opaque = Z_NULL;
			strm.avail_in = 0;
			strm.next_in = Z_NULL;
			//ret = inflateInit(&strm);
			ret = inflateInit2(&strm, 16+MAX_WBITS);
			if (ret != Z_OK) {
				std::cerr << "decompress stream: unkown error" << std::endl;
				exit(1);
			}
			_resetable = false;
			fastBundle.compressType = fc_DECOMPRESSOR;
		}

		void setCompressedFastBundle(FastBundle* pCompressedFastBundle);


		FastBundle* getCompressedFastBundle() {
			if(!_compressedFastBundle) {
				std::cerr << "ERROR " << __FILE__ << " " << __LINE__ << std::endl;
				return NULL;
			}
			else {
				FastBundle* s = _compressedFastBundle;
				_compressedFastBundle = NULL;
				//std::cout << "getCompressedFastBundle" << std::endl;
				return s;
			}

		}

		bool decompress(uint tId);


	};

	/*
	 * extracts reads from each FastBundle and stores them in ReadBundles
	 */
	class FastParser {
	private:
		//////////////////////
		bool _skipEstimate;
		double rerror = 0.0;                   //bella's read error
		////////////////////
		TSeqType _seqType;                          // sequence type
		TFileType _fileType;                        // file type

		uint64 _readsNumber;                        // number of reads
		SyncSwapQueueMPMC<ReadBundle> _syncQueue;   // SyncSwapQueue for ReadBundles

		uint32_t _threadsNumber;                       // number of threads
		std::thread **_processThreads;              // list of threads

		SyncSwapQueueSPSC<FastBundle> **_fastSyncSwapQueues;  // SyncSwapQueue for FastBundles

		FastBundle **_curFastBundles;                         // current FastBundles
		Decompressor **_decompressors;                         // decompressors for each thread

		inline void skipLineBreak(char *&bp, char *&bp_end, const size_t &tId);

		inline void skipLine(char *&bp, char *&bp_end, const size_t &tId);

		inline void skipLine(char *&bp, char *&bp_end, const size_t &l, const size_t &tId);
		
		///////////////////////////
		inline void errorEstimation(char *&bp, char *&bp_end, size_t &l, ReadBundle *&readBundle, ReadBundle *&rbs, const size_t &tId, const char &skip);
		///////////////////////////
		inline void storeLine(
				char *&bp, char *&bp_end, size_t &l,
				ReadBundle *&readBundle, ReadBundle *&rbs, const size_t &tId, const char &skip
		);

		inline void storeSequence(
				char *&bp, char *&bp_end, size_t &l,
				ReadBundle *&readBundle, ReadBundle *&rbs, const size_t &tId, const char &skip
		);

		void nextPart(char *&bp, char *&bp_end, const size_t &tId);

		void processFastq(const size_t &tId);

		void processFasta(const size_t &tId);

		void processMultiline(const size_t &tId);

		StopWatch *_sw;
	public:
		double erate; //for bella's error rate	
		
		double getErate(){
			return this->erate;
		}
	
		SyncSwapQueueMPMC<ReadBundle> *getSyncQueue();          // returns SyncSwapQueue of ReadBundles

		inline uint64 getReadsNumber() { return _readsNumber; } // returns total number of reads

		/*
		 * constructor
		 */
		FastParser(
				uint32 &readBundlesNumber, TFileType fileType, TSeqType seqType,
				SyncSwapQueueSPSC<FastBundle> **_fastSyncSwapQueues,
				const uint32_t &_readerParserThreadsNumber,
				bool skipEstimante
		);

		/*
		 * starts the entire working process
		 */
		void process();

		/*
		 * joins all threads
		 */
		void join();

		void joinWithoutDelete();
		void deleteProcessThread();

		/*
		 * prints some statistical outputs
		 */
		void print();

		/*
		 * destructor
		 */
		~FastParser();
	};

}


#endif /* FASTPARSER_H_ */
