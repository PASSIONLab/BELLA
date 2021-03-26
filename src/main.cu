#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <array>
#include <typeinfo>
#include <tuple>
#include <queue>
#include <memory>
#include <stack>
#include <functional>
#include <cstring>
#include <string.h>
#include <math.h>
#include <cassert>
#include <ios>
#include <chrono>
#include <thread>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <map>
#include <unordered_map>
#include <omp.h>

#include "../include/cxxopts.hpp"

#include "../libcuckoo/cuckoohash_map.hh"
#include "../include/kmercount.hpp"
#include "../include/chain.hpp"
#include "../include/common/bellaio.h"
#include "../include/minimizer.hpp"
#include "../include/syncmer.hpp"

#include "../kmercode/hash_funcs.h"
#include "../kmercode/Kmer.hpp"
#include "../kmercode/Buffer.h"
#include "../kmercode/common.h"
#include "../kmercode/fq_reader.h"
#include "../kmercode/ParallelFASTQ.h"
#include "../kmercode/bound.hpp"

#include "../include/common/utility.h"
#include "../include/common/CSC.h"
#include "../include/common/CSR.h"
#include "../include/common/common.h"
#include "../include/common/IO.h"
#include "../include/overlap.hpp"
#include "../include/align.hpp"

#define LSIZE 16000
#define ITERS 10

#define KMERINDEX uint32_t		
// #define KMERINDEX uint64_t 	// Uncomment to for large genomes and comment out line 60
  
using namespace std;

int main (int argc, char *argv[]) {
	//
	// Program name and purpose
	//
	cxxopts::Options options("BELLA", "Long Read to Long Read Aligner and Overlapper");

	//
	// Setup the input files
	//
	options.add_options()
	("f, fastq", "List of Fastq(s) (required)", 	cxxopts::value<std::string>())
	("o, output", "Output Filename (required)", 	cxxopts::value<std::string>())
	("k, kmer", "K-mer Length", 	 	cxxopts::value<int>()->default_value("17"))
	("x, xdrop", "SeqAn X-Drop", 		cxxopts::value<int>()->default_value("7"))
	("e, error", "Error Rate", 			cxxopts::value<double>()->default_value("0.15"))
	("estimate", "Estimate Error Rate from Data", 			cxxopts::value<bool>()->default_value("false"))
	("skip-alignment", "Overlap Only", 	cxxopts::value<bool>()->default_value("false"))
	("m, memory", "Total RAM of the System in MB", 			cxxopts::value<int>()->default_value("8000"))
	("score-deviation", "Deviation from the Mean Alignment Score [0,1]", 	cxxopts::value<double>()->default_value("0.1"))
	("b, bin-size", "Bin Size for Binning Algorithm", 		cxxopts::value<int>()->default_value("500"))
	("paf", "Output in PAF format", 	cxxopts::value<bool>()->default_value("false"))
	("g, gpus", "GPUs Available", 		cxxopts::value<int>()->default_value("1")) // this must work only if compiled with bella-gpu
	("split-count", "K-mer Counting Split Count", 			cxxopts::value<int>()->default_value("1"))
	("hopc", "Use HOPC representation", cxxopts::value<bool>()->default_value("false"))
	("w, window", "Window Size for Minimizer Selection", 	cxxopts::value<int>()->default_value("0"))
	("s, syncmer", "Enable Syncmer Selection", 				cxxopts::value<bool>()->default_value("false"))
	("u, upper-freq", "K-mer Frequency Upper Bound", 		cxxopts::value<int>()->default_value("8"))
	("l, lower-freq", "K-mer Frequency Lower Bound", 		cxxopts::value<int>()->default_value("2"))
	("h, help", "Usage")
	;

	auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

	char *inputfofn = NULL;	
	if(result.count("fastq")) inputfofn = strdup(result["fastq"].as<std::string>().c_str());
	else
	{
      std::cout << options.help() << std::endl;
      exit(0);		
	}

	char *OutputFile = NULL;	
	if(result.count("output"))
	{
		char* line1 = strdup(result["fastq"].as<std::string>().c_str());
		char* line2 = strdup(".out");

		unsigned int len1 = strlen(line1);
		unsigned int len2 = strlen(line2);

		OutputFile = (char*)malloc(len1 + len2 + 1);
		if (!OutputFile) abort();

		memcpy(OutputFile,  line1, len1);
		memcpy(OutputFile + len1,  line2, len2);
		OutputFile[len1 + len2] = '\0';

		delete line1, line2;
		
		remove(OutputFile);
	}
	else
	{
      std::cout << options.help() << std::endl;
      exit(0);		
	}

	BELLApars bpars;	

	bpars.kmerSize 	= result["kmer"].as<int>();
	bpars.xDrop 	= result["xdrop"].as<int>();
	bpars.errorRate = result["error"].as<double>();

	bpars.estimateErr 	= result["estimate"].as<bool>();
	bpars.skipAlignment = result["skip-alignment"].as<bool>();
	bpars.totalMemory 	= result["memory"].as<int>();

	#ifdef LINUX || OSX // in include/overlap.hpp
		bpars.userDefMem = false;
	#endif

	bpars.deltaChernoff = result["score-deviation"].as<double>();
	if(bpars.deltaChernoff > 1.0 || bpars.deltaChernoff < 0.0)
	{
      std::cout << options.help() << std::endl;
      exit(0);		
	}

	bpars.binSize 	 = result["bin-size"].as<int>();
	bpars.outputPaf	 = result["paf"].as<bool>();
	bpars.numGPU 	 = result["gpus"].as<int>();
	bpars.SplitCount = result["split-count"].as<int>();
	bpars.useHOPC	 = result["hopc"].as<bool>();

	bpars.windowLen  = result["window"].as<int>();
	if(bpars.windowLen != 0)
		bpars.useMinimizer = true;

	bpars.useSyncmer = result["syncmer"].as<bool>();
	if(bpars.useSyncmer)
		bpars.useMinimizer = false;

	int reliableUpperBound	= result["upper-freq"].as<int>();	
	int reliableLowerBound	= result["lower-freq"].as<int>(); 

	// ================ //
	//   Declarations   //
	// ================ //

	vector<filedata> allfiles = GetFiles(inputfofn);
	std::string all_inputs_gerbil = std::string(inputfofn); 
	double ratiophi;
	Kmer::set_k(bpars.kmerSize);
	unsigned int upperlimit = 10000000; // in bytes
	Kmers kmervect;
	vector<string> seqs;
	vector<string> quals;
	vector<string> nametags;
	readVector_ reads;
	Kmers kmersfromreads;

	// vector<tuple<unsigned int, unsigned int, unsigned short int>> occurrences;	// 32 bit, 32 bit, 16 bit (read, kmer, position)
	vector<tuple<KMERINDEX, KMERINDEX, unsigned short int>> transtuples;	// 32 bit, 32 bit, 16 bit (kmer, read, position)

	// ================== //
	// Parameters Summary //
	// ================== //
    
#ifdef PRINT
    printLog(OutputFile);

    std::string kmerSize = std::to_string(bpars.kmerSize);
    printLog(kmerSize);

    std::string GPUs = "ENABLED";
    printLog(GPUs);

	std::string UserDefinedMemory = std::to_string(bpars.totalMemory) + " MB";
	printLog(UserDefinedMemory);

    std::string OutputPAF = std::to_string(bpars.outputPaf);
    printLog(OutputPAF);

    std::string BinSize = std::to_string(bpars.binSize);
    printLog(BinSize);
    
    std::string DeltaChernoff = std::to_string(bpars.deltaChernoff);
    printLog(DeltaChernoff);

    std::string RunPairwiseAlignment = std::to_string(!bpars.skipAlignment);
    printLog(RunPairwiseAlignment);

	if(bpars.useHOPC)
	{
		std::string HOPC = "ENABLED";
		printLog(HOPC);
	}
	else 
	{
		std::string HOPC = "DISABLED";
		printLog(HOPC);
	}

    std::string xDrop = std::to_string(bpars.xDrop);
    printLog(xDrop);

    std::string KmerSplitCount = std::to_string(bpars.SplitCount);
    printLog(KmerSplitCount);

	if(bpars.useMinimizer)
	{
		std::string useMinimizer = "ENABLED";
		printLog(useMinimizer);
		std::string minimizerWindow = std::to_string(bpars.windowLen);
    	printLog(minimizerWindow);
	}
 	else if(bpars.useSyncmer)
	{
		std::string useSyncmer = "ENABLED";
		printLog(useSyncmer);
		std::string minimizerWindow = "0";
    	printLog(minimizerWindow);
	}
	else
	{
		std::string useKmer = "ENABLED";
		printLog(useKmer);
		std::string minimizerWindow = "0";
    	printLog(minimizerWindow);
	}

#endif

	double all;

	// ================ //
	//  Reliable Bound  //
	// ================ //

	int numThreads = 1;
	#pragma omp parallel
	{
	    numThreads = omp_get_num_threads();
	}

	printLog(numThreads);

	// GG: reads global
	vector<readVector_> allreads(MAXTHREADS);

	all = omp_get_wtime();

	// ================ //
	//  K-mer Counting  //
	// ================ //

	CuckooDict<KMERINDEX> countsreliable;

	if(bpars.useSyncmer)
	{
		SyncmerCount(allfiles, countsreliable, reliableLowerBound, reliableUpperBound,
    		upperlimit, bpars);
	}
	else if(bpars.useMinimizer)
    {
    	MinimizerCount(allfiles, countsreliable, reliableLowerBound, reliableUpperBound,
            upperlimit, bpars);
    }
	else
	{
    	SplitCount(allfiles, countsreliable, reliableLowerBound, reliableUpperBound,
               upperlimit, bpars);
	}

	double errorRate;

	if(bpars.useHOPC)
	{
		errorRate = bpars.HOPCerate;
	}
	else
	{
		errorRate = bpars.errorRate;
	}

	printLog(errorRate);
	
	printLog(reliableLowerBound);
	printLog(reliableUpperBound);

	if(bpars.fixedThreshold == -1)
	{
	    ratiophi = slope(bpars.errorRate);
	    float AdaptiveThresholdConstant = ratiophi * (1 - bpars.deltaChernoff);
		printLog(AdaptiveThresholdConstant); 
	}

	// ================ //
	// Fastq(s) Parsing //
	// ================ //

	double parsefastq = omp_get_wtime();

	// vector<vector<tuple<unsigned int, unsigned int, unsigned short int>>> alloccurrences(MAXTHREADS);
	vector<vector<tuple<KMERINDEX, KMERINDEX, unsigned short int>>> alltranstuples(MAXTHREADS);

	unsigned int numReads = 0; // numReads needs to be global (not just per file)

	for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++)
	{
		ParallelFASTQ *pfq = new ParallelFASTQ();
		pfq->open(itr->filename, false, itr->filesize);

		unsigned int fillstatus = 1;
		while(fillstatus)
		{
			fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
			unsigned int nreads = seqs.size();

		#pragma omp parallel for
			for(int i=0; i<nreads; i++) 
			{
				// remember that the last valid position is length()-1
				int len = seqs[i].length();

				readType_ temp;
				nametags[i].erase(nametags[i].begin());	// removing "@"
				temp.nametag = nametags[i];
				temp.seq = seqs[i];    					// save reads for seeded alignment
				temp.readid = numReads+i;
				allreads[MYTHREAD].push_back(temp);
                
                if(bpars.useMinimizer)
                {
                    vector<Kmer> seqkmers;
                    std::vector<int> seqminimizers;    // <position_in_read>
                    for(int j = 0; j <= len - bpars.kmerSize; j++)   // AB: optimize this sliding-window parsing ala HipMer
                    {
                        std::string kmerstrfromfastq = seqs[i].substr(j, bpars.kmerSize);
                        Kmer mykmer(kmerstrfromfastq.c_str(), kmerstrfromfastq.length());
                        seqkmers.emplace_back(mykmer);
                    }

                    getMinimizers(bpars.windowLen, seqkmers, seqminimizers);

                    for(auto minpos: seqminimizers)
                    {
                        std::string strminkmer = seqs[i].substr(minpos, bpars.kmerSize);
                        Kmer myminkmer(strminkmer.c_str(), strminkmer.length());
                        
                        KMERINDEX idx; // kmer_id
                        auto found = countsreliable.find(myminkmer,idx);
                        if(found)
                        {
                            alltranstuples[MYTHREAD].emplace_back(std::make_tuple(idx, numReads+i, minpos));
                        }
                    }
                }
                else
                {
                    for(int j = 0; j <= len - bpars.kmerSize; j++)
                    {
                        std::string kmerstrfromfastq = seqs[i].substr(j, bpars.kmerSize);
                        Kmer mykmer(kmerstrfromfastq.c_str(), kmerstrfromfastq.length());
                        // remember to use only ::rep() when building kmerdict as well
                        Kmer lexsmall;
                        if (bpars.useHOPC)
                        {
                            lexsmall = mykmer.hopc();
                        }
                        else
                        {
                            // remember to use only ::rep() when building kmerdict as well
                            lexsmall = mykmer.rep();
                        }

                        KMERINDEX idx; // kmer_id
                        auto found = countsreliable.find(lexsmall,idx);
                        if(found)
                        {
                            //alloccurrences[MYTHREAD].emplace_back(std::make_tuple(numReads+i, idx, j)); // vector<tuple<numReads,kmer_id,kmerpos>>
                            alltranstuples[MYTHREAD].emplace_back(std::make_tuple(idx, numReads+i, j)); // transtuples.push_back(col_id,row_id,kmerpos)
                        }
                    }
                }
			} // for(int i=0; i<nreads; i++)
			numReads += nreads;
		} //while(fillstatus) 
		delete pfq;

	} // for all files


	KMERINDEX readcount = 0;
	KMERINDEX tuplecount = 0;

	for(int t=0; t<MAXTHREADS; ++t)
	{
		readcount  += allreads[t].size();
		tuplecount += alltranstuples[t].size();
	}

// #define WRITEDATAMATRIX
#ifdef WRITEDATAMATRIX
    WriteToDisk(alltranstuples, countsreliable, readcount, tuplecount);
#endif

	reads.resize(readcount);
	//occurrences.resize(tuplecount);
	transtuples.resize(tuplecount);

	unsigned int readssofar = 0;
	unsigned int tuplesofar = 0;

	for(int t=0; t<MAXTHREADS; ++t)
	{
		copy(allreads[t].begin(), allreads[t].end(), reads.begin()+readssofar);
		readssofar += allreads[t].size();

		//copy(alloccurrences[t].begin(), alloccurrences[t].end(), occurrences.begin() + tuplesofar);
		copy(alltranstuples[t].begin(), alltranstuples[t].end(), transtuples.begin() + tuplesofar);
		tuplesofar += alltranstuples[t].size();
	}

	std::sort(reads.begin(), reads.end());	// bool operator in global.h: sort by readid

	std::vector<string>().swap(seqs);		// free memory of seqs  
	std::vector<string>().swap(quals);		// free memory of quals

	std::string fastqParsingTime = std::to_string(omp_get_wtime() - parsefastq) + " seconds";
	printLog(fastqParsingTime);
	printLog(numReads);

	// ====================== //
	// Sparse Matrix Creation //
	// ====================== //
	
	unsigned int nkmer = countsreliable.size();
	
	// to help the parsing script
    cout << nkmer << endl;
	
	double matcreat = omp_get_wtime();
	CSC<KMERINDEX, unsigned short int> transpmat(transtuples, nkmer, numReads,
							[] (unsigned short int& p1, unsigned short int& p2) 
							{
								return p1;
							}, false);	// hashspgemm doesn't require sorted rowids within each column
	// remove memory of transtuples
	std::vector<tuple<KMERINDEX, KMERINDEX, unsigned short int>>().swap(transtuples);

	std::string TransposeSparseMatrixCreationTime = std::to_string(omp_get_wtime() - matcreat) + " seconds";
	printLog(TransposeSparseMatrixCreationTime);


	double transbeg = omp_get_wtime();	
	CSC<KMERINDEX, unsigned short int> spmat = transpmat.Transpose();
	std::string ReTransposeTime = std::to_string(omp_get_wtime() - transbeg) + " seconds";
	printLog(ReTransposeTime);


	// ==================================================== //
	// Sparse Matrix Multiplication (aka Overlap Detection) //
	// ==================================================== //
		
	spmatPtr_ getvaluetype(make_shared<spmatType_>());
	HashSpGEMMGPU(
		spmat, transpmat, 
		// n-th k-mer positions on read i and on read j
	    [&bpars, &reads] (const unsigned short int& begpH, const unsigned short int& begpV, 
	        const unsigned int& id1, const unsigned int& id2)
		{
			spmatPtr_ value(make_shared<spmatType_>());

			std::string& read1 = reads[id1].seq;
			std::string& read2 = reads[id2].seq;

			// GG: function in chain.h
			multiop(value, read1, read2, begpH, begpV, bpars.kmerSize);
			return value;
		},
	    [&bpars, &reads] (spmatPtr_& m1, spmatPtr_& m2, const unsigned int& id1, 
	        const unsigned int& id2)
		{
			// GG: after testing correctness, these variables can be removed
			std::string& readname1 = reads[id1].nametag;
			std::string& readname2 = reads[id2].nametag;

			// GG: function in chain.h
			chainop(m1, m2, bpars, readname1, readname2);
			return m1;
		},
	    reads, getvaluetype, OutputFile, bpars, ratiophi);

	double totaltime = omp_get_wtime()-all;

    std::string TotalRuntime = std::to_string(totaltime) + " seconds";   
    printLog(TotalRuntime);
	
	cout << totaltime << endl;

	return 0;
}
