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

#ifndef KMCWRITER_H_
#define KMCWRITER_H_

#include "SyncQueue.h"
#include "Bundle.h"

namespace gerbil {

class KmcWriter {
	std::string _fileName;
	std::thread* _processThread;
	FILE* _file;
	uint32_t _k;
	TOutputFormat _outputFormat;
	std::vector<std::pair<std::string,uint32> > *listKmer;

	SyncSwapQueueMPSC<KmcBundle>* _kmcSyncSwapQueue;

	uint64_t _fileSize;
public:
	KmcWriter(std::string fileName,
			SyncSwapQueueMPSC<KmcBundle>* kmcSyncSwapQueue, const uint32_t &k, const TOutputFormat pOutputFormat);
	~KmcWriter();
	
	std::vector<std::pair<std::string,uint32>> *getListKmer(){
		return this->listKmer;
	}

	void process();
	void join();
	void joinWithoutDelete();
	void deleteProcessThread();
	void print();
};

}

#endif /* KMCWRITER_H_ */
