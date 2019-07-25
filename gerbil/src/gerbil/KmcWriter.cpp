/*
 * KmcWriter.cpp
 *
 *  Created on: 07.06.2015
 *      Author: marius
 */

#include "../../include/gerbil/KmcWriter.h"
#include<tuple>
#include<vector>

gerbil::KmcWriter::KmcWriter(std::string fileName, SyncSwapQueueMPSC<KmcBundle>* kmcSyncSwapQueue, const uint32_t &k, const TOutputFormat pOutputFormat) {
	_processThread = NULL;
	_fileName = fileName;
	_k = k;
	_outputFormat = pOutputFormat;
	_kmcSyncSwapQueue = kmcSyncSwapQueue;
	listKmer = new std::vector<std::pair<std::string,uint32> >();
	if(_outputFormat != of_none) {
		std::remove(_fileName.c_str());
		_file = fopen(_fileName.c_str(), "wb");
		if (!_file) {
			std::cerr << "unable to create output-file" << std::endl;
			exit(21);
		}
	}
	if(_outputFormat == of_gerbil)
		setbuf(_file, NULL);
	_fileSize = 0;
}

gerbil::KmcWriter::~KmcWriter() {

}


void gerbil::KmcWriter::process() {
	
	if(_processThread)
		return;
	_processThread = new std::thread([this]{
		IF_MESS_KMCWRITER(
			StopWatch sw;
			sw.start();
		)
		IF_DEB(
			printf("kmcWriter start...\n");
		)
		KmcBundle* kb = new KmcBundle;
		std::string kmer;
		if(_outputFormat == of_fasta) {
			uint32 counter;
			char kmerSeq[_k + 1]; kmerSeq[_k] = '\0';
			const size_t kMerSize_B = (_k + 3) / 4;
			const char c[4] = {'A', 'C', 'G', 'T'};

			IF_MESS_KMCWRITER(sw.hold();)
			while(_kmcSyncSwapQueue->swapPop(kb)) {
				IF_MESS_KMCWRITER(sw.proceed();)
				if(!kb->isEmpty()) {
					const byte* p = kb->getData();
					const byte* end = p + kb->getSize();
					while(p < end) {
						counter = (uint32)*(p++);
						if(counter >= 255) {
							// large value
							counter = *((uint32*)p);
							p += 4;
						}
						// k-mer: convert bytes to string
						for(uint i = 0; i < _k; ++i)
							kmerSeq[i] = c[(p[i >> 2] >> (2 * (3 - (i & 0x3)))) & 0x3];

						// increase pointer
						p += kMerSize_B;
						kmer = kmerSeq;
						
						listKmer->push_back(std::make_pair(kmer,counter));
												
		// print fasta (console/file)
					//	fprintf(_file, ">%u\n%s\n", counter, kmerSeq);
					}
				}
				kb->clear();
				IF_MESS_KMCWRITER(sw.hold();)
			}
			_fileSize = ftell(_file);
		}
		else if(_outputFormat == of_gerbil) {
			IF_MESS_KMCWRITER(sw.hold();)
			while(_kmcSyncSwapQueue->swapPop(kb)) {
				_fileSize += kb->getSize();
				IF_MESS_KMCWRITER(sw.proceed();)
				if(!kb->isEmpty()) {
					fwrite ((char*) kb->getData() , 1 , kb->getSize() , _file );
				}
				kb->clear();
				IF_MESS_KMCWRITER(sw.hold();)
			}
		}
		IF_MESS_KMCWRITER(sw.proceed();)

		delete kb;
		if(_file)
			fclose(_file);
		IF_MESS_KMCWRITER(
			sw.stop();
			printf("kmcWriter: %7.3f s\n", sw.get_s());
		)
	});
}

void gerbil::KmcWriter::join() {
	_processThread->join();
	delete _processThread;
}

void gerbil::KmcWriter::joinWithoutDelete(){
	_processThread->join();
}

void gerbil::KmcWriter::deleteProcessThread(){
	delete _processThread;
}


void gerbil::KmcWriter::print() {
	printf("size of output  : %12.3f MB\n", (double)_fileSize / 1024 / 1024);
}
