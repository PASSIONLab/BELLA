/*
 * Application.cpp
 *
 *  Created on: 10.06.2015
 *      Author: marius
 */

#include "../../include/gerbil/Application.h"
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#ifdef __linux__
#include <sys/sysinfo.h>
#elif _WIN32
#include <windows.h>
#else
// is there another OS?
#endif

#ifdef GPU
#include <cuda_runtime.h>
#endif


/**
 * Detect the total amount of memory at the system.
 * @return
 */
unsigned long long getTotalSystemMemory()  {
#ifdef __linux__
	struct sysinfo info;
	if(sysinfo(&info) == -1)
		return DEF_MEMORY_SIZE;
	return info.totalram;
#elif _WIN32
	MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys;
#else
	return DEF_MEMORY_SIZE;
#endif
}


/**
 * Detect the amount of free memory at the system.
 * @return
 */
unsigned long long getFreeSystemMemory()  {
#ifdef __linux__
	struct sysinfo info;
	if(sysinfo(&info) == -1)
		return DEF_MEMORY_SIZE;
	return info.freeram;
#elif _WIN32
	MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullAvailPhys;
#else
	return DEF_MEMORY_SIZE;
#endif
}



/*
 * Constructor
 * default for params
 */
// TODO: Reorder attributes
gerbil::Application::Application(uint32_t kmerSize, 
				std::string fastFileName,
				std::string tempFolderName, 
				uint32_t thresholdMin,
				std::string kmcFileName,
				bool skipEstimate) :
		_k(kmerSize), _m(0), _tempFilesNumber(0), _sequenceSplitterThreadsNumber(0),
		_superSplitterThreadsNumber(0), _hasherThreadsNumber(0), _thresholdMin(thresholdMin), _memSize(0),
		_threadsNumber(0), _norm(DEF_NORM),
		_fastFileName(fastFileName), _tempFolderName(tempFolderName), _kmcFileName(kmcFileName), _tempFiles(NULL),
		_rtRun1(0.0), _rtRun2(0.0), _memoryUsage1(0), _memoryUsage2(0),
		_readerParserThreadsNumber(1), _numGPUs(0),
		_singleStep(0), _leaveBinStat(false), _histogram(0), _outputFormat(of_fasta), _skipEstimate(skipEstimate)
{

}

gerbil::Application::~Application() {
#ifdef GPU
	cudaDeviceReset();
#endif
}


void gerbil::Application::process() {
	verbose = true;
	autocompleteParams();
	
	printParamsInfo();
	
	checkParams();	

	ReadBundle::setK(_k);

	if (_singleStep != 2)
		run1();
	else
		loadBinStat();

	if (_singleStep != 1)
		run2();

	if(_singleStep != 1 && !_leaveBinStat)
		std::remove((_tempFolderName + "binStatFile.txt").c_str());
}

void gerbil::Application::run1() {
	// no buffer for output stream
	setbuf(stdout, NULL);

	// time
	StopWatch sw(CLOCK_REALTIME);
	sw.start();

	// calculate memory
	uint32 frBlocksNumber, readBundlesNumber, superBundlesNumber;
	uint64 superWriterBufferSize;
	distributeMemory1(frBlocksNumber, readBundlesNumber, superBundlesNumber,
			superWriterBufferSize);



	// init pipeline
	FastReader fastReader(frBlocksNumber, _fastFileName,
			_readerParserThreadsNumber);
	FastParser fastParser(readBundlesNumber, fastReader.getFileType(), st_reads,
			fastReader.getSyncSwapQueues(), _readerParserThreadsNumber, _skipEstimate);
	SequenceSplitter sequenceSplitter(superBundlesNumber,
			fastParser.getSyncQueue(), _sequenceSplitterThreadsNumber, _k, _m,
			_tempFilesNumber, _norm);
	SuperWriter superWriter(_tempFolderName,
			sequenceSplitter.getSuperBundleQueues(), _tempFilesNumber,
			superWriterBufferSize);

	// start pipepline
	fastReader.process();
	fastParser.process();
	sequenceSplitter.process();
	superWriter.process();

	// join all
	fastReader.join();
	//fastParser.join();
	//joins all the threads, retrive the value and then delete the threads
	fastParser.joinWithoutDelete();
	erate = fastParser.getErate();
	fastParser.deleteProcessThread();
	
	//std::cout<<erate<<'\n';
	sequenceSplitter.join();
	superWriter.join();

	// save bin statistic
	_tempFiles = superWriter.getTempFiles();
	saveBinStat();

	sw.stop();
	_rtRun1 = sw.get_s();
	//return error rate for bella

	// verbose output
	if (verbose) {
		printf("================== STAGE 1 ==================\n");
		fastReader.print();
		fastParser.print();
		sequenceSplitter.print();
		superWriter.print();
		printf("memory usage           : %12lu MB\n", B_TO_MB(_memoryUsage1));
		printf("---------------------------------------------\n");
	}
}

void gerbil::Application::run2() {
	// no buffer for output stream
	setbuf(stdout, NULL);

	// time
	StopWatch sw(CLOCK_REALTIME);
	sw.start();

	// calculate statistic for tempFiles
	TempFileStatistic tempFileStatistic(_tempFiles, _tempFilesNumber);

	// distribute number of hashers to number of gpu and cpu hashers
	uint8_t hasherThreadsGPU = std::min(_numGPUs, _hasherThreadsNumber);
	uint8_t hasherThreadsCPU = std::max(0,_hasherThreadsNumber-hasherThreadsGPU);

	// calculate memory
	uint32 superBundlesNumber, kmcBundlesNumber, kMerBundlesNumber;
	uint64 maxKmcHashtableSize;
	distributeMemory2(&tempFileStatistic, superBundlesNumber, kmcBundlesNumber,
			maxKmcHashtableSize, kMerBundlesNumber);

	// init distributor
	KmerDistributer distributor(hasherThreadsCPU, hasherThreadsGPU, _tempFilesNumber);

	// init pipeline
	SuperReader superReader(superBundlesNumber, _tempFiles, _tempFilesNumber, &distributor);
	KmerHasher kmerHasher(_k, kmcBundlesNumber, superReader.getSuperBundleQueue(),
			_superSplitterThreadsNumber, hasherThreadsCPU, hasherThreadsGPU,
			_tempFiles, _tempFilesNumber, _thresholdMin, _norm, _tempFolderName,
			maxKmcHashtableSize, kMerBundlesNumber,
			superReader.getTempFilesOrder(), &distributor);
	KmcWriter kmcWriter(_kmcFileName, kmerHasher.getKmcSyncSwapQueue(), _k, _outputFormat);

	// start pipeline
	superReader.process();
	kmerHasher.process();
	kmcWriter.process();

	// join all
	superReader.join();
	kmerHasher.join();
	//kmcWriter.join();

	//join, retrive the value and then delete the process
	kmcWriter.joinWithoutDelete();
	listKmer = kmcWriter.getListKmer();
	kmcWriter.deleteProcessThread();
	

	sw.stop();
	_rtRun2 = sw.get_s();

	// save bin statistic
	saveBinStat();

	if(_histogram)
		kmerHasher.saveHistogram();
		
	// verbose output
	if (verbose) {
		printf("================== STAGE 2 ==================\n");
		kmerHasher.print();
		kmcWriter.print();
		printf("memory usage    : %12lu MB\n", B_TO_MB(_memoryUsage2));
		printf("---------------------------------------------\n");
		printSummary();
	}
}

/*
 * private
 */

void gerbil::Application::distributeMemory1(uint32 &fastBundlesNumber,
		uint32 &readBundlesNumber, uint32 &superBundlesNumber,
		uint64 &superWriterBufferSize) {
	// some space for general consumption
	uint64 base_memory_B = RUN1_MEMORY_GENERAL_B;

	// const memory for SequenceSplitter
	base_memory_B += (_tempFilesNumber * _sequenceSplitterThreadsNumber)
			* SUPER_BUNDLE_DATA_SIZE_B;

	_memoryUsage1 = base_memory_B;

	// FastReader
	fastBundlesNumber = MIN_FASTBUNDLEBUFFER_SIZE_B / FAST_BUNDLE_DATA_SIZE_B;
	_memoryUsage1 += fastBundlesNumber * FAST_BUNDLE_DATA_SIZE_B;

	// FastParser
	readBundlesNumber = MIN_READBUNDLEBUFFER_SIZE_B / READ_BUNDLE_SIZE_B;
	_memoryUsage1 += readBundlesNumber * READ_BUNDLE_SIZE_B;

	// SequenceSplitter
	superBundlesNumber =
	MIN_SUPERBUNDLEBUFFER_SIZE_B / SUPER_BUNDLE_DATA_SIZE_B;
	_memoryUsage1 += superBundlesNumber * SUPER_BUNDLE_DATA_SIZE_B;

	// SuperWriter
	_memoryUsage1 += MIN_SUPERWRITERBUFFER_SIZE_B;

	// check minimal memory
	if (_memoryUsage1 > MB_TO_B(_memSize)) {
		std::cerr
				<< "Not enough memory!\n\t Reduced number of temporary files/ the Number of SequenceSplitters or increasing the amount of memory."
				<< std::endl;
		exit(1);
	}

	// compute memory for buffers

	// memory for FastBundles
	uint64 optFastBundlesNumber = MAX_FASTBUNDLEBUFFER_SIZE_B
			/ FAST_BUNDLE_DATA_SIZE_B;
	uint64 memOptFastBundles = 0;
	if (optFastBundlesNumber > fastBundlesNumber) {
		optFastBundlesNumber -= fastBundlesNumber; // already assured
		memOptFastBundles = optFastBundlesNumber * FAST_BUNDLE_DATA_SIZE_B;
	} else
		optFastBundlesNumber = 0;

	// memory for ReadBundles
	uint64 optReadBundlesNumber = MAX_READBUNDLEBUFFER_SIZE_B
			/ READ_BUNDLE_SIZE_B;
	uint64 memOptReadBundles = 0;
	if (optReadBundlesNumber > readBundlesNumber) {
		optReadBundlesNumber -= readBundlesNumber; // already assured
		memOptReadBundles = optReadBundlesNumber * READ_BUNDLE_SIZE_B;
	} else
		optReadBundlesNumber = 0;

	// memory for SuperBundles
	uint64 optSuperBundlesNumber = MAX_SUPERBUNDLEBUFFER_SIZE_B
			/ SUPER_BUNDLE_DATA_SIZE_B;
	uint64 memOptSuperBundles = 0;
	if (optSuperBundlesNumber > superBundlesNumber) {
		optSuperBundlesNumber -= superBundlesNumber; // already assured
		memOptSuperBundles = optSuperBundlesNumber * SUPER_BUNDLE_DATA_SIZE_B;
	} else
		optSuperBundlesNumber = 0;

	uint64 availableMemory_B = (MB_TO_B(_memSize) - _memoryUsage1) * 0.5; // assure 50% memory for SuperWriterBuffer

	// assure memory
	uint64 sumOptMem = memOptFastBundles + memOptReadBundles
			+ optSuperBundlesNumber;
	if (sumOptMem <= availableMemory_B) {
		fastBundlesNumber += optFastBundlesNumber;
		readBundlesNumber += optReadBundlesNumber;
		superBundlesNumber += optSuperBundlesNumber;
	} else {
		fastBundlesNumber += (double) availableMemory_B
				/ FAST_BUNDLE_DATA_SIZE_B
				* ((double) memOptFastBundles / sumOptMem);
		readBundlesNumber += (double) availableMemory_B / READ_BUNDLE_SIZE_B
				* ((double) memOptReadBundles / sumOptMem);
		superBundlesNumber += (double) availableMemory_B
				/ SUPER_BUNDLE_DATA_SIZE_B
				* ((double) memOptSuperBundles / sumOptMem);
	}

	_memoryUsage1 = base_memory_B;
	_memoryUsage1 += fastBundlesNumber * FAST_BUNDLE_DATA_SIZE_B;
	_memoryUsage1 += readBundlesNumber * READ_BUNDLE_SIZE_B;
	_memoryUsage1 += superBundlesNumber * SUPER_BUNDLE_DATA_SIZE_B;

	availableMemory_B = MB_TO_B(_memSize) - _memoryUsage1;

	superWriterBufferSize = availableMemory_B
			/ (SUPER_BUNDLE_DATA_SIZE_B + sizeof(SuperBundleStackItem));

	_memoryUsage1 += superWriterBufferSize
			* (SUPER_BUNDLE_DATA_SIZE_B + sizeof(SuperBundleStackItem));
}

void gerbil::Application::distributeMemory2(
		TempFileStatistic* tempFileStatistic, uint32 &superBundlesNumber,
		uint32 &kmcBundlesNumber, uint64 &maxKmcHashtableSize,
		uint32 &kMerBundlesNumber) {
	uint64 bytesPerHashEntry = getKMerByteNumbers(_k) + sizeof(uint_cv);

	// some space for general consumption
	uint64 base_memory_B = RUN2_MEMORY_GENERAL_B;

	// const memory consumption
	base_memory_B += _numGPUs * GPU_COPY_BUFFER_SIZE;

	base_memory_B += _superSplitterThreadsNumber * _hasherThreadsNumber
			* KMER_BUNDLE_DATA_SIZE_B;
	base_memory_B += 2 * (_hasherThreadsNumber + _numGPUs) * KMER_BUNDLE_DATA_SIZE_B;
	base_memory_B += 2 * (_hasherThreadsNumber + _numGPUs)
			* (1 + FAILUREBUFFER_KMER_BUNDLES_NUMBER_PER_THREAD)
			* KMER_BUNDLE_DATA_SIZE_B;
	base_memory_B += (_hasherThreadsNumber + _numGPUs) * KMC_BUNDLE_DATA_SIZE_B;

	_memoryUsage2 = base_memory_B;
	uint64 availableMemory_B = MB_TO_B(_memSize) - _memoryUsage2;

	// compute minimum sizes
	superBundlesNumber = MIN_SUPERBUNDLEBUFFER2_SIZE_B
			/ SUPER_BUNDLE_DATA_SIZE_B;
	kMerBundlesNumber = MIN_KMERBUNDLEBUFFER_SIZE_B / KMER_BUNDLE_DATA_SIZE_B;
	kmcBundlesNumber = MIN_KMCBUNDLEBUFFER_SIZE_B / KMC_BUNDLE_DATA_SIZE_B;
	maxKmcHashtableSize = MIN_KMCHASHTABLE_SIZE_B / bytesPerHashEntry;

	// update memory usage for step 2
	_memoryUsage2 += (uint64) superBundlesNumber * SUPER_BUNDLE_DATA_SIZE_B;
	_memoryUsage2 += (uint64) kMerBundlesNumber * KMER_BUNDLE_DATA_SIZE_B;
	_memoryUsage2 += (uint64) kmcBundlesNumber * KMC_BUNDLE_DATA_SIZE_B;
	_memoryUsage2 += maxKmcHashtableSize * bytesPerHashEntry;
	availableMemory_B = MB_TO_B(_memSize) - _memoryUsage2;
	//printf("%u\t\t%u\t\t%u\t\t%lu\n", superBundlesNumber, kMerBundlesNumber, kmcBundlesNumber, maxKmcHashtableSize);

	// compute memory for hashtable
	uint64 maxUKMersNumber = std::min(tempFileStatistic->getMaxKMersNumber(),
			tempFileStatistic->getAvg2SdKMersNumber());
	if (maxUKMersNumber > maxKmcHashtableSize) {
		maxUKMersNumber -= maxKmcHashtableSize; // already assured
		uint64 extraSize;
		if (maxUKMersNumber
				* bytesPerHashEntry<= availableMemory_B * MEM_KEY_HT)
			extraSize = maxUKMersNumber;
		else
			extraSize = availableMemory_B * MEM_KEY_HT / bytesPerHashEntry;
		maxKmcHashtableSize += extraSize;
		_memoryUsage2 += extraSize * bytesPerHashEntry;
		availableMemory_B = MB_TO_B(_memSize) - _memoryUsage2;
	}

	// compute memory for buffers
	uint64 optSuperBundlesNumber =
			tempFileStatistic->getAvg2SdSize() / SUPER_BUNDLE_DATA_SIZE_B;
	uint64 memOptSuperBundles = 0;
	if (optSuperBundlesNumber > superBundlesNumber) {
		optSuperBundlesNumber -= superBundlesNumber; // already assured
		memOptSuperBundles = optSuperBundlesNumber * SUPER_BUNDLE_DATA_SIZE_B;
	} else
		optSuperBundlesNumber = 0;

	uint64 optKMerBundlesNumber = tempFileStatistic->getAvg2SdKMersNumber()
			* getKMerByteNumbers(_k) / KMER_BUNDLE_DATA_SIZE_B;
	//printf("optkmer: %lu\n", optKMerBundlesNumber * KMER_BUNDLE_DATA_SIZE_B);
	uint64 memOptKMerBundles = 0;
	if (optKMerBundlesNumber > kMerBundlesNumber) {
		optKMerBundlesNumber -= kMerBundlesNumber; // already assured
		memOptKMerBundles = optKMerBundlesNumber * KMER_BUNDLE_DATA_SIZE_B;
	} else
		optKMerBundlesNumber = 0;

	uint64 optKmcBundlesNumber = tempFileStatistic->getAvgKMersNumber()
			/ (1 + 2 * log((uint64) _thresholdMin))
			* getKMerByteNumbers(_k) / KMC_BUNDLE_DATA_SIZE_B; // should be: number of uk-Mers
			//printf("optkmc: %lu\n", optKmcBundlesNumber* KMC_BUNDLE_DATA_SIZE_B);
	uint64 memOptKmcBundles = 0;
	if (optKmcBundlesNumber > kmcBundlesNumber) {
		optKmcBundlesNumber -= kmcBundlesNumber; // already assured
		memOptKmcBundles = optKmcBundlesNumber * KMC_BUNDLE_DATA_SIZE_B;
	} else
		optKmcBundlesNumber = 0;

	// assure memory
	uint64 sumOptMem = memOptSuperBundles + memOptKMerBundles
			+ memOptKmcBundles;
	if (sumOptMem <= availableMemory_B) {
		superBundlesNumber += optSuperBundlesNumber;
		kMerBundlesNumber += optKMerBundlesNumber;
		kmcBundlesNumber += optKmcBundlesNumber;
	} else {
		superBundlesNumber += (double) availableMemory_B
				/ SUPER_BUNDLE_DATA_SIZE_B
				* ((double) memOptSuperBundles / sumOptMem);
		kMerBundlesNumber += (double) availableMemory_B
				/ KMER_BUNDLE_DATA_SIZE_B
				* ((double) memOptKMerBundles / sumOptMem);
		kmcBundlesNumber += (double) availableMemory_B / KMC_BUNDLE_DATA_SIZE_B
				* ((double) memOptKmcBundles / sumOptMem);
	}

	// update memory usage for step 2
	_memoryUsage2 = base_memory_B;
	_memoryUsage2 += superBundlesNumber * SUPER_BUNDLE_DATA_SIZE_B;
	_memoryUsage2 += kMerBundlesNumber * KMER_BUNDLE_DATA_SIZE_B;
	_memoryUsage2 += kmcBundlesNumber * KMC_BUNDLE_DATA_SIZE_B;
	_memoryUsage2 += maxKmcHashtableSize * bytesPerHashEntry;
	availableMemory_B = MB_TO_B(_memSize) - _memoryUsage2;
	//printf("%u\t\t%u\t\t%u\t\t%lu\n", superBundlesNumber, kMerBundlesNumber, kmcBundlesNumber, maxKmcHashtableSize);
	/*printf("%lu MB\t\t%lu MB\t\t%lu MB\t\t%lu MB\n",
	 B_TO_MB(superBundlesNumber * SUPER_BUNDLE_DATA_SIZE_B),
	 B_TO_MB(kMerBundlesNumber * KMER_BUNDLE_DATA_SIZE_B),
	 B_TO_MB(kmcBundlesNumber * KMC_BUNDLE_DATA_SIZE_B),
	 B_TO_MB(maxKmcHashtableSize * bytesPerHashEntry));*/

	//printf("avg Size     = %12lu\n", tempFileStatistic->getAvgKMersNumber());
	//printf("stabw        = %12lu\n", tempFileStatistic->getSdKMersNumber());
	//printf("max Size     = %12lu\n", tempFileStatistic->getMaxKMersNumber());
	//printf("avg+2Sd Size = %12lu\n", tempFileStatistic->getAvg2SdKMersNumber());
}

void gerbil::Application::checkSystem() {
	printf("check System...\n");
	printf("supports lockfree\n");
	printf("\tuint64 : ");
	if (__atomic_always_lock_free(sizeof(uint64), 0))
		printf("yes\n");
	else
		printf("no\n");
	printf("\tuint32 : ");
	if (__atomic_always_lock_free(sizeof(uint32), 0))
		printf("yes\n");
	else
		printf("no\n");
	printf("\tuint16 : ");
	if (__atomic_always_lock_free(sizeof(uint16), 0))
		printf("yes\n");
	else
		printf("no\n");
	printf("\tuint8  : ");
	if (__atomic_always_lock_free(sizeof(uint8), 0))
		printf("yes\n");
	else
		printf("no\n");
	unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
	printf("supported threads: %u\n", concurentThreadsSupported);
	//int sysinfo(struct sysinfo *info);
	/*printf("total-ram: %lu\n", (size_t)sysconf( _SC_PHYS_PAGES ) *
	 (size_t)sysconf( _SC_PAGESIZE ));*/

	//struct sysinfo info;
	//sysinfo(&info);
	uint64_t totalram = getTotalSystemMemory();
	uint64_t freeram = getFreeSystemMemory();
	printf("total-ram: %2.3f GB\n",
			(double) totalram / 1024 / 1024 / 1024);
	printf("free-ram: %2.3f GB\n", (double) freeram / 1024 / 1024 / 1024);

	/*struct stat64 stat_buf;
	if (!_fastFileName.empty() && !stat64(_fastFileName.c_str(), &stat_buf)) {
		printf("file size: %f GB\n",
				(double) stat_buf.st_size / 1024 / 1024 / 1024);

	}*/
}

void gerbil::Application::autocompleteParams() {
#define SET_DEFAULT(x, d) if(!x) x = d

	// set to default values
	SET_DEFAULT(_k, DEF_KMER_SIZE);
	SET_DEFAULT(_threadsNumber,
			std::thread::hardware_concurrency() < MIN_THREADS_NUMBER ? DEF_THREADS_NUMBER : std::thread::hardware_concurrency());
	if (_threadsNumber < 4)
		_threadsNumber = 4;
	SET_DEFAULT(_sequenceSplitterThreadsNumber,
			_threadsNumber <= 4 ? 2 : _threadsNumber - 3);

#if false
	SET_DEFAULT(_superSplitterThreadsNumber, 1 + (_threadsNumber - 1) * 4 / 10);
	SET_DEFAULT(_hasherThreadsNumber, 1 + (_threadsNumber - 1) * 8 / 10);
#endif

	SET_DEFAULT(_superSplitterThreadsNumber, 1 + (_threadsNumber - 1) * 5 / 10);
	SET_DEFAULT(_hasherThreadsNumber, _numGPUs + 1 + (_threadsNumber - 1) * 6 / 10);
	
	SET_DEFAULT(_tempFilesNumber, DEF_TEMPFILES_NUMBER);
	SET_DEFAULT(_thresholdMin, DEF_THRESHOLD_MIN);

	// set size of available memory
	if (!_memSize) {
		//struct sysinfo info;
		//if (!sysinfo(&info)) {
		uint64_t totalram = getTotalSystemMemory();
			_memSize = totalram / 1024 / 1024;
			if (_memSize > 1024)				// for system and file buffers
				_memSize -= 1024;
			else
				_memSize = 0;					// no auto setup
		//} else
		//	_memSize = DEF_MEMORY_SIZE;
	}

	if (!_sequenceSplitterThreadsNumber)
		_sequenceSplitterThreadsNumber = 1;
	if (!_superSplitterThreadsNumber)
		_superSplitterThreadsNumber = 1;
	if (!_hasherThreadsNumber)
		_hasherThreadsNumber = 1;

	// set size of minimizers
	if (!_m) {
		uint_tfn x = _tempFilesNumber;
		while (x) {
			x >>= 2;
			++_m;
		}
		// for balanced distribution of file sizes
		if (_m < MAX_MINIMIZER_SIZE)
			++_m;
		if (_m < MAX_MINIMIZER_SIZE)
			++_m;

		if (_m < MIN_MINIMIZER_SIZE)
			_m = MIN_MINIMIZER_SIZE;

	}
}

void gerbil::Application::checkParams() {
#define CHECK_PARAM_IN_RANGE(x, min, max, s1, s2)																			\
		if (x < min) {																											\
			std::cerr << s1 << " (" << (uint64)x << ") " << " is too small (should be: " << s2 << " >= " << min << ")" << std::endl;		\
			exit(1);																											\
		}																														\
		else if (x > max) {																										\
			std::cerr << s1 << " (" << (uint64)x << ") " << " is too large (should be: " << s2 << " <= " << max << ")" << std::endl;		\
			exit(1);																											\
		}

	// check file paths
	if (_fastFileName.empty() && _singleStep != 2) {
		printf("input file is unknown\n");
		exit(1);
	}
	if (_tempFolderName.empty()) {
		printf("temp folder is unknown\n");
		exit(1);
	}
	if (_kmcFileName.empty() && _singleStep != 1) {
		printf("output file is unknown\n");
		exit(1);
	}

	// check parameter ranges
	CHECK_PARAM_IN_RANGE(_memSize, MIN_MEMORY_SIZE, MAX_MEMORY_SIZE,
			"size e of memory", "e");
	CHECK_PARAM_IN_RANGE(_tempFilesNumber, MIN_TEMPFILES_NUMBER,
			MAX_TEMPFILES_NUMBER, "number n of temp files", "n");
	CHECK_PARAM_IN_RANGE(_m, MIN_MINIMIZER_SIZE, MAX_MINIMIZER_SIZE,
			"size m of minimizers", "m");
	CHECK_PARAM_IN_RANGE(_k, MIN_KMER_SIZE, MAX_KMER_SIZE, "size k of k-mers",
			"k");
	CHECK_PARAM_IN_RANGE(_thresholdMin, MIN_THRESHOLD_MIN, MAX_THRESHOLD_MIN,
			"minimum counts l of k-mer", "l");

	CHECK_PARAM_IN_RANGE(_threadsNumber, MIN_THREADS_NUMBER, MAX_THREADS_NUMBER,
			"number t of threads", "t");

	// check additional conditions
	if (_m > _k) {
		std::cerr
				<< "size m of minimizers is too large, size of k-mers k is too small (should be: m <= k)\n";
		exit(1);
	}
	if ((((uint64) 4) << (2 * _m)) < _tempFilesNumber) {
		std::cerr
				<< "number n of temp files is too large, size of minimizers m is too small (should be: n <= 4^m)\n";
		exit(1);
	}
}

void gerbil::Application::printParamsInfo() {
	printf("______________________________________________\n");
	printf("Gerbil version %i.%i\n", VERSION_MAJOR, VERSION_MINOR);
	printf("================= PARAMETERS =================\n");
	printf("size of k-mers          : %5u\n", _k);
	printf("size of minimizers      : %5u\n", _m);
	printf("threshold min           : %5u\n", _thresholdMin);
	printf("normalized kmers        : ");
	_norm ? printf(" true\n") : printf("false\n");
	printf("number of temp-files    : %5lu\n", _tempFilesNumber);
	printf("total number of threads : %5u\n", _threadsNumber);
	printf("number of splitters     : %5u\n", _sequenceSplitterThreadsNumber);
	printf("number of hashers       : %5u\n", _hasherThreadsNumber);
	printf("input                   :       %s\n", _fastFileName.c_str());
	printf("temp                    :       %s\n", _tempFolderName.c_str());
	printf("output                  :       %s\n", _kmcFileName.c_str());
	printf("size of memory          : %5lu MB\n", _memSize);
	printf("number of gpu's         : %5u\n", _numGPUs);
	printf("---------------------------------------------\n");
}

void gerbil::Application::printSummary() {
	printf("===================================================\n");
	printf("================== S U M M A R Y ==================\n");
	printf("===================================================\n");
	printf("time (real) for stage1    :  % 5.3f s\n", _rtRun1);
	printf("time (real) for stage2    :  % 5.3f s\n", _rtRun2);
	printf("time (real) total (1 + 2) :  % 5.3f s\n", _rtRun1 + _rtRun2);
	printf("---------------------------------------------------\n");
}

void gerbil::Application::saveBinStat() {
	std::string binStatFileName(_tempFolderName);
	binStatFileName += "binStatFile.txt";
	FILE* binStatFile = fopen(binStatFileName.c_str(), "wb+");
	for (uint_tfn tempFileId(0); tempFileId < _tempFilesNumber; ++tempFileId)
		_tempFiles[tempFileId].fprintStat(binStatFile);
	fclose(binStatFile);
}

void gerbil::Application::loadBinStat() {
	if (_tempFiles)
		return;
	std::string binStatFileName(_tempFolderName);
	binStatFileName += "binStatFile.txt";
	FILE* binStatFile = fopen(binStatFileName.c_str(), "rb");
	_tempFiles = new TempFile[_tempFilesNumber];
	for (uint_tfn tempFileId(0); tempFileId < _tempFilesNumber; ++tempFileId) {
		_tempFiles[tempFileId].loadStats(_tempFolderName, binStatFile);
	}
	fclose(binStatFile);
}
