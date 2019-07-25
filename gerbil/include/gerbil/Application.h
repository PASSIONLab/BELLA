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

#ifndef APPLICATION_H_
#define APPLICATION_H_

#include <vector>

#include <string>

#include "global.h"
#include "FastReader.h"
#include "FastParser.h"
#include "SequenceSplitter.h"
#include "SuperWriter.h"
#include "SuperReader.h"
#include "KmerHasher.h"
#include "KmcWriter.h"
#include "TempFileStatistic.h"

namespace gerbil {

/*
 * Prog
 */
class Application {
	double erate;                                                //error rate
	std::vector<std::pair<std::string,uint32>> *listKmer;	     //list of the counted kmer
	bool _skipEstimate;                                           //execute or skip the error estimation
	uint32_t _k;							// size of k-mers
	uint8 _m;								// size of minimizers
	uint_tfn _tempFilesNumber;				// number of tempFiles
	uint8 _sequenceSplitterThreadsNumber;	// number of threads for SequenceSplitter
	uint8 _superSplitterThreadsNumber;		// number of threads for SuperSplitter
	uint8 _hasherThreadsNumber;				// number of hash/extract threads
	std::string _fastFileName;				// filename of fast[a/q] (with path)
	std::string _tempFolderName;			// foldername of temp (with path)
	std::string _kmcFileName;				// filename of kmc (with path)
	uint32 _thresholdMin;					// min k-mer counter to store
	uint64 _memSize;						// size of ram in MB
	uint8 _threadsNumber;					// total number of threads
	uint8 _readerParserThreadsNumber;		// number of readers and parsers
	uint8_t _numGPUs;						// number of gpu's to use
	TOutputFormat  _outputFormat;           // which output format to use
	bool _norm;								// normalization of kmers enabled
	uint _singleStep;						// processes only one step (default: 0 => all steps)
	bool _leaveBinStat;						// leaves binStatFile
	bool _histogram;						// prints histogram

	TempFile* _tempFiles;					// bin files (step1 --> step2)

	double _rtRun1;							// realtime for run1
	double _rtRun2;							// realtime for run2

	uint64 _memoryUsage1;					// memory usage for step1
	uint64 _memoryUsage2;					// memory usage for step2

	void checkSystem();

	void autocompleteParams();
	void checkParams();
	void printParamsInfo();

	void printSummary();

	void distributeMemory1(uint32 &frBlocksNumber, uint32 &readBundlesNumber,
			uint32 &superBundlesNumber, uint64 &superWriterBufferSize);
	void distributeMemory2(TempFileStatistic* tempFileStatistic,
			uint32 &superBundlesNumber, uint32 &kmcBundlesNumber,
			uint64 &maxKmcHashtableSize, uint32 &kMerBundlesNumber);

	void saveBinStat();
	void loadBinStat();
public:
	Application(uint32_t kmerSize,
			std::string fastFileName,
			std::string tempFolderName,
			uint32_t thresholdMin,
			std::string kmcFileName,
			bool skipEstimate);
	~Application();

	double getErate(){
		return this->erate;
	}

	std::vector<std::pair<std::string,unsigned int>> *getListKmer(){
		return this->listKmer;
	}
	


	void process();
	void run1();

	void run2();
};

}

#endif /* APPLICATION_H_ */
