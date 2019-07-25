/*
 * FastParser.cpp
 *
 *  Created on: 21.05.2015
 *      Author: marius
 */

#include <iostream>
#include <fstream>

#include "../../include/gerbil/FastParser.h"
#define ACSCIIBASE 33


void gerbil::Decompressor::setCompressedFastBundle(FastBundle *pCompressedFastBundle) {
	if (_compressedFastBundle)
		std::cerr << "ERROR" << __FILE__ << " " << __LINE__ << std::endl;
	else {
		_compressedFastBundle = pCompressedFastBundle;
		strm.avail_in = _compressedFastBundle->size;
		strm.next_in = (byte *) _compressedFastBundle->data;
	}
}


bool gerbil::Decompressor::decompress(uint tId) {
	if (!strm.avail_in)
        return _compressedFastBundle->size != 0;
	strm.avail_out = FAST_BUNDLE_DATA_SIZE_B - fastBundle.size;
	strm.next_out = ((byte *) fastBundle.data) + fastBundle.size;
	ret = inflate(&strm, Z_NO_FLUSH);
	// ignore
	assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
	switch (ret) {
        case Z_OK: break;
        case Z_STREAM_END:
            if(strm.avail_in != 0) {
                fprintf(stderr, "ERROR: unexpected end of stream (thread[%u])! Archive corrupt?!\n", tId);
                exit(66);
            }
            break;
        default:
            (void) inflateEnd(&strm);
            fprintf(stderr, "ERROR: unhandled return value ret[%u] = %i (Archive corrupt?)\n", tId, ret);
            exit(66);
	}

	fastBundle.size = FAST_BUNDLE_DATA_SIZE_B - strm.avail_out;
	return fastBundle.size < FAST_BUNDLE_DATA_SIZE_B;
}

void gerbil::FastParser::nextPart(char *&bp, char *&bp_end, const size_t &tId) {
	FastBundle *curFastBundle = _curFastBundles[tId];
	bool finish = false;
	if (curFastBundle->compressType == fc_none || curFastBundle->size == 0) {
		curFastBundle->clear();
		finish = !_fastSyncSwapQueues[tId]->swapPop(curFastBundle);
		_curFastBundles[tId] = curFastBundle;
		if (!finish && curFastBundle->compressType != fc_none) {
			_decompressors[tId]->reset();
			_decompressors[tId]->setCompressedFastBundle(curFastBundle);
		}
	}

	if (finish) {
		return;
	}

	if (curFastBundle->compressType == fc_none) {
		bp = curFastBundle->data;
		bp_end = curFastBundle->data + curFastBundle->size;
	} else {
		if (curFastBundle->compressType == fc_DECOMPRESSOR) {
			_decompressors[tId]->fastBundle.clear();
		}
		while (_decompressors[tId]->decompress(tId)) {
			// aufraeumen
			curFastBundle = _decompressors[tId]->getCompressedFastBundle();
			curFastBundle->clear();
			if (_fastSyncSwapQueues[tId]->swapPop(curFastBundle)) {
				_decompressors[tId]->setCompressedFastBundle(curFastBundle);
			} else {
				std::cerr << "ERROR " << __FILE__ << "  " << __LINE__ << std::endl;
				break;
			}
		}

		if (curFastBundle->size == 0) {
			_curFastBundles[tId] = _decompressors[tId]->getCompressedFastBundle();
			_curFastBundles[tId]->compressType = fc_none;
			_decompressors[tId]->fastBundle.finalize(fc_DECOMPRESSOR);
		} else
			_curFastBundles[tId] = &(_decompressors[tId]->fastBundle);
		bp = _decompressors[tId]->fastBundle.data;
		bp_end = bp + _decompressors[tId]->fastBundle.size;
    }
}


inline void gerbil::FastParser::skipLineBreak(char *&bp, char *&bp_end, const size_t &tId) {
	if (!*(++bp))
		nextPart(bp, bp_end, tId);
	if (*bp == '\n' && !*(++bp))
		nextPart(bp, bp_end, tId);
}

inline void gerbil::FastParser::skipLine(char *&bp, char *&bp_end, const size_t &tId) {
	while (*bp != '\n' && *bp != '\r')
		if (!*(++bp))
			nextPart(bp, bp_end, tId);
	skipLineBreak(bp, bp_end, tId);
}

inline void gerbil::FastParser::skipLine(char *&bp, char *&bp_end, const size_t &l, const size_t &tId) {
	if (bp + l < bp_end)
		bp += l; //skip +\n
	else {
		size_t skip = l - (bp_end - bp);
		nextPart(bp, bp_end, tId);
		bp += skip;
	}
	skipLineBreak(bp, bp_end, tId);
}
////////////////////////////////////////
inline void gerbil::FastParser::errorEstimation(char *&bp, char *&bp_end, size_t &l, ReadBundle *&readBundle, ReadBundle *&rbs, const size_t &tId, const char &skip){
	rerror = 0.0;
	int length = 0;
	int bqual;
	double berror;
	do{
		char *sl(bp);
		//error rate calculation
		while(*bp != '\n' && *bp != '\r' && *bp){
			length++;	
			bqual = (int)*bp - ACSCIIBASE;
			berror = pow(10, -(double)bqual/10);
			rerror += berror;
			++bp;
		}
		//part from the storeSequence 
		if(!*bp){
			nextPart(bp, bp_end, tId);
		}
		if(*bp && (*bp== '\n' || *bp == '\r')){
			++bp;
			if(!*bp){
			nextPart(bp, bp_end, tId);
			}
		}

	} while(*bp && skip && *bp != skip);
	rerror = rerror/length;
	
	

}
////////////////////////////////////////



inline void gerbil::FastParser::storeSequence(
		char *&bp, char *&bp_end, size_t &l,
		ReadBundle *&readBundle, ReadBundle *&rbs, const size_t &tId, const char &skip
) {
	bool first = true;
	do {
		char *sl(bp);
		while (*bp != '\n' && *bp != '\r' && *bp) {
			++bp;
		}
		l = bp - sl;

		// store first part
		if (first) {
			if (!readBundle->add(l, sl)) {
				//readBundle->print();
				_syncQueue.swapPush(readBundle);
				readBundle->add(l, sl);
			}
		} else {
			if (!readBundle->expand(l, sl)) {
				readBundle->transferKm1(rbs);
				//readBundle->print();
				_syncQueue.swapPush(readBundle);
				ReadBundle *swap = rbs;
				rbs = readBundle;
				readBundle = swap;
				readBundle->expand(l, sl);
			}
		}
		if (!*bp)
			nextPart(bp, bp_end, tId);
		if(*bp && (*bp == '\n' || *bp == '\r')) {
			++bp;
			if (!*bp) {
				nextPart(bp, bp_end, tId);
			}
		}
		first = false;
	} while (*bp && skip && *bp != skip);
}


inline void gerbil::FastParser::storeLine(
		char *&bp, char *&bp_end, size_t &l,
		ReadBundle *&readBundle, ReadBundle *&rbs, const size_t &tId, const char &skip
) {

	bool first = true;

	do {

		char *sl(bp);
		while (*bp != '\n' && *bp != '\r' && *bp) {
			++bp;
		}
		l = bp - sl;

		// store first part
		if (first) {
			if (!readBundle->add(l, sl)) {
				_syncQueue.swapPush(readBundle);
				readBundle->add(l, sl);
			}
		} else {
			if (!readBundle->expand(l, sl)) {
				readBundle->transfer(rbs);
				_syncQueue.swapPush(readBundle);
				ReadBundle *swap = rbs;
				rbs = readBundle;
				readBundle = swap;
				readBundle->expand(l, sl);
			}
		}

		if (!*bp) {
			nextPart(bp, bp_end, tId);
			// store second part
			sl = bp;
			while (*bp != '\n' && *bp != '\r') ++bp;
			l += bp - sl;
			if (!readBundle->expand(bp - sl, sl)) {
				readBundle->transfer(rbs);
				_syncQueue.swapPush(readBundle);
				ReadBundle *swap = rbs;
				rbs = readBundle;
				readBundle = swap;
				readBundle->expand(bp - sl, sl);
			}
		}
		//readBundle->print();
		skipLineBreak(bp, bp_end, tId);
		first = false;
	} while (*bp && skip && *bp != skip);
}



gerbil::FastParser::FastParser(
		uint32 &readBundlesNumber, TFileType fileType,
		TSeqType seqType, SyncSwapQueueSPSC<FastBundle> **fastSyncSwapQueues,
		const uint32_t &_readerParserThreadsNumber,
		bool skipEstimate
) : _syncQueue(readBundlesNumber), _threadsNumber(_readerParserThreadsNumber), _processThreads(NULL) {
	_fileType = fileType;
	_seqType = seqType;
	_fastSyncSwapQueues = fastSyncSwapQueues;
	///////
	_skipEstimate = skipEstimate;
	//////
	_readsNumber = 0;

	_curFastBundles = new FastBundle *[_threadsNumber];
	_decompressors = new Decompressor *[_threadsNumber];
	_processThreads = new std::thread *[_threadsNumber];
	IF_MESS_FASTPARSER(
			_sw = new StopWatch[_threadsNumber];
			for (uint32_t i(0); i < _threadsNumber; ++i)
				_sw[i].setMode(CLOCK_THREAD_CPUTIME_ID);
	)
}

void gerbil::FastParser::processFastq(const size_t &tId) {
	char *bp;
	char *bp_end;
	size_t l;
	
	ReadBundle *readBundle = new ReadBundle();
	ReadBundle *rbs = new ReadBundle();
	//////
	erate = 0.0;
	//////
	while (true) {
		nextPart(bp, bp_end, tId);
		if (!_curFastBundles[tId]->size)
			break;
		while (*bp) {
			// skip description
			skipLine(bp, bp_end, tId);

			// store read
			storeSequence(bp, bp_end, l, readBundle, rbs, tId, '+');

			// skip + [description]
			skipLine(bp, bp_end, tId);

			// skip qualifiers
			if(_skipEstimate){
            			skipLine(bp, bp_end, l, tId);
			}else{
				errorEstimation(bp, bp_end,l,readBundle, rbs, tId,'@');
			}

			erate+=rerror;
			while (*bp && *bp != '@')
                		skipLine(bp, bp_end, tId);

			++_readsNumber;
		}
	}
	erate =(double) erate/_readsNumber;
        // std::cout<<erate<<'\n';		

	
	if (!readBundle->isEmpty())
		_syncQueue.swapPush(readBundle);
	delete readBundle;
	delete rbs;
}

void gerbil::FastParser::processFasta(const size_t &tId) {
	char *bp;
	char *bp_end;
	size_t l;

	ReadBundle *readBundle = new ReadBundle();
	ReadBundle *rbs = new ReadBundle();

	while (true) {
		nextPart(bp, bp_end, tId);
		while (*bp) {
			// skip description
			skipLine(bp, bp_end, tId);

			storeSequence(bp, bp_end, l, readBundle, rbs, tId, '>');

			++_readsNumber;
		}
		if (!_curFastBundles[tId]->size && !*bp)
			break;
	}
	if (!readBundle->isEmpty()) {
		_syncQueue.swapPush(readBundle);
	}
	delete readBundle;
	delete rbs;
}

void gerbil::FastParser::processMultiline(const size_t &tId) {
	char *bp;
	char *bp_end;
	size_t l;

	ReadBundle *readBundle = new ReadBundle();
	ReadBundle *rbs = new ReadBundle();

	while (true) {
		nextPart(bp, bp_end, tId);
		if (!_curFastBundles[tId]->size)
			break;
		while (*bp) {

			storeLine(bp, bp_end, l, readBundle, rbs, tId, 0);

			++_readsNumber;
		}
	}
	if (!readBundle->isEmpty())
		_syncQueue.swapPush(readBundle);
	delete readBundle;
	delete rbs;
}


void gerbil::FastParser::process() {
	if (_seqType != st_reads) {
		std::cerr << "unsupported sequence type for fastq" << std::endl;
		exit(1);
	}
	for (uint32_t i = 0; i < _threadsNumber; i++) {
		_processThreads[i] = new std::thread([this](uint64_t tId) {
			IF_MESS_FASTPARSER(_sw[tId].start();)

			_curFastBundles[tId] = new FastBundle();
			_decompressors[tId] = new Decompressor();

			switch (_fileType) {
				case ft_fastq:
					processFastq(tId);
					break;
				case ft_fasta:
					processFasta(tId);
					break;
				case ft_multiline:
					processMultiline(tId);
					break;
				default:
					std::cerr << "unknown Filetype" << std::endl;
					exit(1);
			}

			delete _curFastBundles[tId];
			delete _decompressors[tId];

			IF_MESS_FASTPARSER(
					_sw[tId].stop();
					printf("time parser[%2lu]: %.3f s\n", tId, _sw[tId].get_s());
			)
		}, i);
	}
}

void gerbil::FastParser::join() {
	for (uint32_t i = 0; i < _threadsNumber; ++i) {
		_processThreads[i]->join();
		delete _processThreads[i];
	}
	_syncQueue.finalize();
	//printf("fastParser is rdy...\n");
}

void gerbil::FastParser::joinWithoutDelete(){
	for(uint32 i = 0; i<_threadsNumber; ++i){
		_processThreads[i]->join();
	}
}

void gerbil::FastParser::deleteProcessThread(){
	for(uint32 i = 0; i<_threadsNumber; ++i){
		delete _processThreads[i];
	}
	_syncQueue.finalize();
}

gerbil::SyncSwapQueueMPMC<gerbil::ReadBundle> *gerbil::FastParser::getSyncQueue() {
	return &_syncQueue;
}

void gerbil::FastParser::print() {
	printf("number of reads        : %12lu\n", _readsNumber);
}

gerbil::FastParser::~FastParser() {
	delete[] _curFastBundles;
	delete[] _decompressors;
	delete[] _processThreads;
	IF_MESS_FASTPARSER(delete[] _sw;)
}

