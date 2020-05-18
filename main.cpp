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

#include "libcuckoo/cuckoohash_map.hh"
#include "kmercount.h"
#include "chain.h"
#include "bellaio.h"

#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include "kmercode/Buffer.h"
#include "kmercode/common.h"
#include "kmercode/fq_reader.h"
#include "kmercode/ParallelFASTQ.h"
#include "kmercode/bound.hpp"

#include "mtspgemm2017/utility.h"
#include "mtspgemm2017/CSC.h"
#include "mtspgemm2017/CSR.h"
#include "mtspgemm2017/common.h"
#include "mtspgemm2017/IO.h"
#include "mtspgemm2017/overlapping.h"
#include "mtspgemm2017/align.h"

#define LSIZE 16000
#define ITERS 10

#define KMERINDEX uint32_t		
// #define KMERINDEX uint64_t 	// Uncomment to for large genomes and comment out line 56
  
using namespace std;

int main (int argc, char *argv[]) {
	//
	// Program name and purpose
	//
	cout << "\nBELLA: Long Read to Long Read Aligner and Overlapper\n" << endl;
	//
	// Setup the input files
	//
	option_t *optList, *thisOpt;
	// Get list of command line options and their arguments 
	// Follow an option with a colon to indicate that it requires an argument.

	optList = NULL;
	optList = GetOptList(argc, argv, (char*)"f:o:c:d:hk:a:ze:x:c:m:r:ps:qg:u:w:l:i:");

	char	*all_inputs_fofn 	= NULL;	// List of fastqs (i)
	char	*OutputFile 		= NULL;	// output filename (o)
	int		InputCoverage 		= 0;	// Coverage required (d)

	BELLApars b_parameters;

	if(optList == NULL)
	{
        std::string ErrorMessage("BELLA execution terminated: not enough parameters or invalid option. Run with -h to print out the command line options.\n");
		printLog(ErrorMessage);
		return 0;
	}

	while (optList!=NULL) {
		thisOpt = optList;
		optList = optList->next;
		switch (thisOpt->option) {
			case 'f': {
				if(thisOpt->argument == NULL)
				{
					std::string ErrorMessage = "BELLA execution terminated: -f requires an argument. Run with -h to print out the command line options.\n";
					printLog(ErrorMessage);
					return 0;
				}
				all_inputs_fofn = strdup(thisOpt->argument);
				break;
			}
			case 'p': b_parameters.outputPaf = true; break; // PAF format
			case 'o': {
				if(thisOpt->argument == NULL)
				{
					std::string ErrorMessage = "BELLA execution terminated: -o requires an argument. Run with -h to print out the command line options.\n";
					printLog(ErrorMessage);
					return 0;
				}

				char* line1 = strdup(thisOpt->argument);
				char* line2 = strdup(".out");

				unsigned int len1 = strlen(line1);
				unsigned int len2 = strlen(line2);

				OutputFile = (char*)malloc(len1 + len2 + 1);
				if (!OutputFile) abort();

				memcpy(OutputFile, line1, len1);
				memcpy(OutputFile + len1, line2, len2);
				OutputFile[len1 + len2] = '\0';

				delete line1;
                delete line2;
                
                // Delete file to avoid errors in output
                remove(OutputFile);

				break;
			}
			case 'c': {
				if(thisOpt->argument == NULL)
				{
					std::string ErrorMessage = "BELLA execution terminated: -c requires an argument. Run with -h to print out the command line options.\n";
					printLog(ErrorMessage);
					return 0;
				}
				InputCoverage = atoi(thisOpt->argument);  
				break;
			}
			case 'z': b_parameters.skipAlignment = true; break;
			case 'k': {
				b_parameters.kmerSize = atoi(thisOpt->argument);
				break;
			}
			case 'r': {
				b_parameters.minProbability = stod(thisOpt->argument);
				break;
			}
			case 'e': { // User suggests error rate
				if(thisOpt->argument == NULL)
				{
					std::string ErrorMessage = "BELLA execution terminated: -e requires an argument. Run with -h to print out the command line options.\n";
					printLog(ErrorMessage);
				}	
				b_parameters.errorRate = strtod(thisOpt->argument, NULL);
				break;
			}
			case 'q': { // Read set has quality values and BELLA can use it to estimate the error rate
				b_parameters.skipEstimate = false;
				break;
			}
			case 'u': {	// Default: skipEstimate and errorRate = 0.15
				b_parameters.errorRate = 0.15;	// Default value
				break;
			}
			case 'a': {
				b_parameters.fixedThreshold = atoi(thisOpt->argument);
				break;
			}
			case 'x': {
				b_parameters.xDrop = atoi(thisOpt->argument);
				break;
			}
			case 'w': {
				b_parameters.binSize = atoi(thisOpt->argument);
				break;
			}
            case 'g': {
				b_parameters.numGPU = atoi(thisOpt->argument);
				break;
			}
			case 'm': {
				b_parameters.totalMemory = stod(thisOpt->argument);
				b_parameters.userDefMem = true;
				std::string UserDefinedMemory = std::to_string(b_parameters.totalMemory) + " MB";
				printLog(UserDefinedMemory);
				break;
			}
			case 's': {
				b_parameters.SplitCount = atoi(thisOpt->argument);
				break;
			}
			case 'd': {
				if(stod(thisOpt->argument) > 1.0 || stod(thisOpt->argument) < 0.0)
				{
					std::string ErrorMessage = "BELLA execution terminated: -d requires a value in [0, 1]. Run with -h to print out the command line options.\n";
					printLog(ErrorMessage);
					return 0;
				}
				b_parameters.deltaChernoff = stod(thisOpt->argument);
				break;
			}
			case 'h': {
				cout << "Usage:\n" << endl;
				cout << "	-f : List of fastq(s)	(required)" 	<< endl;
				cout << "	-o : Output filename	(required)" 	<< endl;
				cout << "	-c : Dataset InputCoverage	(required)" << endl;
				cout << "	-k : KmerSize [17]" 					<< endl;
				cout << "	-a : User-defined alignment threshold [FALSE, -1]" 		<< endl;
				cout << "	-x : SeqAn xDrop [7]" 									<< endl;
				cout << "	-e : Error rate [0.15]" 				<< endl;
				cout << "	-q : Estimare error rate from the dataset [FALSE]" 	<< endl;
				cout << "	-u : Use default error rate setting [FALSE]"		<< endl;
				cout << "	-m : Total RAM of the system in MB [auto estimated if possible or 8,000 if not]"	<< endl;
				cout << "	-z : Do not run pairwise alignment [FALSE]" 			<< endl;
				cout << "	-d : Deviation from the mean alignment score [0.10]"	<< endl;
				cout << "	-w : Bin size binning algorithm [500]" 	<< endl;
				cout << "	-p : Output in PAF format [FALSE]" 		<< endl;
				cout << "	-r : Probability threshold for reliable range [0.002]"  << endl;
                cout << "	-g : GPUs available [1, only works when BELLA is compiled for GPU]" 	<< endl;
				cout << "	-s : K-mer counting split count can be increased for large dataset [1]\n" 	<< endl;

				FreeOptList(thisOpt); // Done with this list, free it
				return 0;
			}
		}
	}

	if(all_inputs_fofn == NULL || OutputFile == NULL || InputCoverage == 0)
	{
		std::string ErrorMessage = "BELLA execution terminated: missing arguments. Run with -h to print out the command line options.\n";
		printLog(ErrorMessage);

		return 0;
    }
    
	if(b_parameters.errorRate == 0.00 && b_parameters.skipEstimate == true)
	{
		std::string str1 = "BELLA execution terminated. The user should either:\n\n";
		std::string str2 = "	* -e = suggest an error rate;\n";
		std::string str3 = "	* -q = confirm that the data has quality values and we can estimate the error rate from the data set;\n";
		std::string str4 = "	* -u = confirm that we can use a default error rate (0.15).\n"; // This might not be worth it for large runs (diBELLA)
		std::string ErrorMessage = str1 + str2 + str3 + str4;

		printLog(ErrorMessage);

		return 0;
	}

	free(optList);
    free(thisOpt);
    
	// ================ //
	// 	 Declarations   //
	// ================ //
    
	vector<filedata> allfiles = GetFiles(all_inputs_fofn);
	std::string all_inputs_gerbil = std::string(all_inputs_fofn); 
	int reliableLowerBound, reliableUpperBound; // reliable range reliableLowerBound and reliableUpperBound bound
	double ratiophi;
	Kmer::set_k(b_parameters.kmerSize);
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
    printLog(InputCoverage);

    std::string kmerSize = std::to_string(b_parameters.kmerSize);
    printLog(kmerSize);

	std::string GPUs = "DISABLED";
    printLog(GPUs);

    std::string OutputPAF = std::to_string(b_parameters.outputPaf);
    printLog(OutputPAF);

    std::string BinSize = std::to_string(b_parameters.binSize);
    printLog(BinSize);
    
    std::string DeltaChernoff = std::to_string(b_parameters.deltaChernoff);
    printLog(DeltaChernoff);

    std::string RunPairwiseAlignment = std::to_string(!b_parameters.skipAlignment);
    printLog(RunPairwiseAlignment);

	if(b_parameters.fixedThreshold == -1)
	{
		std::string AdaptiveAlignmentThreshold = "ENABLED";
		printLog(AdaptiveAlignmentThreshold);
	}
	else 
	{
		std::string AdaptiveAlignmentThreshold = "DISABLE";
		std::string FixedAlignmentThreshold = std::to_string(b_parameters.fixedThreshold);
    	printLog(AdaptiveAlignmentThreshold);
		printLog(FixedAlignmentThreshold);
	}

    std::string xDrop = std::to_string(b_parameters.xDrop);
    printLog(xDrop);

    std::string ReliableCutoffProbability = std::to_string(b_parameters.minProbability);
    printLog(ReliableCutoffProbability);

 	std::string KmerSplitCount = std::to_string(b_parameters.SplitCount);
    printLog(KmerSplitCount);

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

	SplitCount(allfiles, countsreliable, reliableLowerBound, reliableUpperBound, 
		InputCoverage, upperlimit, b_parameters);

	double errorRate  = b_parameters.errorRate;
	printLog(errorRate);
	printLog(reliableLowerBound);
	printLog(reliableUpperBound);

	if(b_parameters.fixedThreshold == -1)
	{
	    ratiophi = slope(b_parameters.errorRate);
	    float AdaptiveThresholdConstant = ratiophi * (1 - b_parameters.deltaChernoff);
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

				for(int j = 0; j <= len - b_parameters.kmerSize; j++)  
				{
					std::string kmerstrfromfastq = seqs[i].substr(j, b_parameters.kmerSize);
					Kmer mykmer(kmerstrfromfastq.c_str(), kmerstrfromfastq.length());
					// remember to use only ::rep() when building kmerdict as well
					Kmer lexsmall = mykmer.rep();

					KMERINDEX idx; // kmer_id
					auto found = countsreliable.find(lexsmall,idx);
					if(found)
					{
						//alloccurrences[MYTHREAD].emplace_back(std::make_tuple(numReads+i, idx, j)); // vector<tuple<numReads,kmer_id,kmerpos>>
						alltranstuples[MYTHREAD].emplace_back(std::make_tuple(idx, numReads+i, j)); // transtuples.push_back(col_id,row_id,kmerpos)
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
	HashSpGEMM(
		spmat, transpmat, 
		// n-th k-mer positions on read i and on read j
	    [&b_parameters, &reads] (const unsigned short int& begpH, const unsigned short int& begpV, 
	        const unsigned int& id1, const unsigned int& id2)
		{
			spmatPtr_ value(make_shared<spmatType_>());

			std::string& read1 = reads[id1].seq;
			std::string& read2 = reads[id2].seq;

			// GG: function in chain.h
			multiop(value, read1, read2, begpH, begpV, b_parameters.kmerSize);
			return value;
		},
	    [&b_parameters, &reads] (spmatPtr_& m1, spmatPtr_& m2, const unsigned int& id1, 
	        const unsigned int& id2)
		{
			// GG: after testing correctness, these variables can be removed
			std::string& readname1 = reads[id1].nametag;
			std::string& readname2 = reads[id2].nametag;

			// GG: function in chain.h
			chainop(m1, m2, b_parameters, readname1, readname2);
			return m1;
		},
	    reads, getvaluetype, OutputFile, b_parameters, ratiophi);

    std::string TotalRuntime = std::to_string(omp_get_wtime()-all) + " seconds";   
    printLog(TotalRuntime);

	return 0;
}
