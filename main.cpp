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
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <map>
#include <unordered_map>
#include <omp.h>

#include "libcuckoo/cuckoohash_map.hh"
#include "kmercount.h"
#include "chain.h"

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
#define PRINT
//#define JELLYFISH

double deltaChernoff = 0.1;            

using namespace std;

int main (int argc, char *argv[]) {

	//
	// Program name and purpose
	//
	cout << "\nBELLA - Long Read Aligner for De Novo Genome Assembly\n" << endl;
	//
	// Setup the input files
	//
	option_t *optList, *thisOpt;
	// Get list of command line options and their arguments 
	// Follow an option with a colon to indicate that it requires an argument.

	optList = NULL;
	optList = GetOptList(argc, argv, (char*)"f:i:o:d:hk:a:ze:x:c:m:r:pb:s:q:gu:");

	char	*kmer_file 			= NULL;	// Reliable k-mer file from Jellyfish
	char	*all_inputs_fofn 	= NULL;	// List of fastqs (i)
	char	*out_file 			= NULL;	// output filename (o)
	int		coverage 			= 0;	// Coverage required (d)

	BELLApars b_parameters;

	if(optList == NULL)
	{
		cout << "BELLA execution terminated: not enough parameters or invalid option" << endl;
		cout << "Run with -h to print out the command line options\n" << endl;
		return 0;
	}

	while (optList!=NULL) {
		thisOpt = optList;
		optList = optList->next;
		switch (thisOpt->option) {
			case 'f': {
#ifdef JELLYFISH
				if(thisOpt->argument == NULL)
				{
					cout << "BELLA execution terminated: -f requires an argument" << endl;
					cout << "Run with -h to print out the command line options\n" << endl;
					return 0;
				}
#endif
				kmer_file = strdup(thisOpt->argument);
				break;
			}
			case 'i': {
				if(thisOpt->argument == NULL)
				{
					cout << "BELLA execution terminated: -i requires an argument" << endl;
					cout << "Run with -h to print out the command line options\n" << endl;
					return 0;
				}
				all_inputs_fofn = strdup(thisOpt->argument);
				break;
			}
			case 'p': b_parameters.outputPaf = true; break; // PAF format
			case 'o': {
				if(thisOpt->argument == NULL)
				{
					cout << "BELLA execution terminated: -o requires an argument" << endl;
					cout << "Run with -h to print out the command line options\n" << endl;
					return 0;
				}
				char* line1 = strdup(thisOpt->argument);
				char* line2 = strdup(".out");
				size_t len1 = strlen(line1);
				size_t len2 = strlen(line2);

				out_file = (char*)malloc(len1 + len2 + 1);
				if (!out_file) abort();

				memcpy(out_file, line1, len1);
				memcpy(out_file + len1, line2, len2);
				out_file[len1 + len2] = '\0';

				delete line1;
				delete line2;

				break;
			}
			case 'd': {
				if(thisOpt->argument == NULL)
				{
					cout << "BELLA execution terminated: -d requires an argument" << endl;
					cout << "Run with -h to print out the command line options\n" << endl;
					return 0;
				}
				coverage = atoi(thisOpt->argument);  
				break;
			}
			case 'z': b_parameters.skipAlignment = true; break;
			case 'k': {
				b_parameters.kmerSize = atoi(thisOpt->argument);
				break;
			}
			case 'r': {
				b_parameters.kmerRift = atoi(thisOpt->argument);
				break;
			}
			case 's': {
				b_parameters.minSurvivedKmers = atoi(thisOpt->argument);
				break;
			}
			case 'e': { // User suggests erro rate
				b_parameters.skipEstimate = true;
				b_parameters.errorRate = strtod(thisOpt->argument, NULL);
				break;
			}
			case 'q': { // Read set has quality values and BELLA can use it to estimate the error rate
				b_parameters.skipEstimate = false;
				break;
			}
			case 'u': {	// Default: skipEstimate and errorRate = 0.15
				b_parameters.skipEstimate = true;
				b_parameters.errorRate = 0.15;	// Default value
				break;
			}
			 case 'g': { // use Gerbil as kmerCounter
				b_parameters.useGerbil = true;
				break;
			}
			case 'a': {
				b_parameters.defaultThr = atoi(thisOpt->argument);
				b_parameters.adapThr = false;
				break;
			}
			case 'x': {
				b_parameters.xDrop = atoi(thisOpt->argument);
				break;
			}
			case 'b': {
				b_parameters.binSize = atoi(thisOpt->argument);
				break;
			}
			case 'm': {
				b_parameters.totalMemory = stod(thisOpt->argument);
				b_parameters.userDefMem = true;
				cout << "User defined memory set to " << b_parameters.totalMemory << " MB " << endl;
				break;
			}
			case 'c': {
				if(stod(thisOpt->argument) > 1.0 || stod(thisOpt->argument) < 0.0)
				{
					cout << "BELLA execution terminated: -c requires a value in [0, 1]" 	<< endl;
					cout << "Run with -h to print out the command line options\n" 		<< endl;
					return 0;
				}
				b_parameters.deltaChernoff = stod(thisOpt->argument);
				break;
			}
			case 'h': {
				cout << "Usage:\n" << endl;
				cout << " -f : List from Jellyfish (required if Jellyfish kmerCounting is used)" 	<< endl; // Reliable k-mers are selected by BELLA
				cout << " -i : List of fastq(s)	(required)" 	<< endl;
				cout << " -o : Output filename	(required)" 	<< endl;
				cout << " -d : Dataset coverage	(required)" 	<< endl; // TO DO: add coverage estimation
				cout << " -k : KmerSize [17]" 					<< endl;
				cout << " -a : User-defined alignment threshold [false, 0]" 		<< endl;
				cout << " -x : SeqAn xDrop [7]" 									<< endl;
				cout << " -e : Error rate [0.15]" 				<< endl;
				cout << " -q : Estimare error rate from the dataset [false]" 	<< endl;
				cout << " -u : Use default error rate setting [false]"			<< endl;
				cout << " -g : Use Gerbil as kmerCounter [false]" 				<< endl;
				cout << " -m : Total RAM of the system in MB [auto estimated if possible or 8,000 if not]" << endl;
				cout << " -z : Do not run pairwise alignment [false]" 				<< endl;
				cout << " -c : Deviation from the mean alignment score [0.10]" 		<< endl;
				cout << " -r : KmerRift: bases separating two k-mers [kmerSize]"	<< endl;
				cout << " -s : Common k-mers threshold to compute alignment [1]"	<< endl;
				cout << " -b : Bin size binning algorithm [500]" 	<< endl;
				cout << " -p : Output in PAF format [false]\n" 		<< endl;

				FreeOptList(thisOpt); // Done with this list, free it
				return 0;
			}
		}
	}

#ifdef JELLYFISH
	if(kmer_file == NULL || all_inputs_fofn == NULL || out_file == NULL || coverage == 0)
	{
		cout << "BELLA execution terminated: missing arguments" << endl;
		cout << "Run with -h to print out the command line options\n" << endl;
		return 0;
	}
#else
	if(all_inputs_fofn == NULL || out_file == NULL || coverage == 0)
	{
		cout << "BELLA execution terminated: missing arguments" << endl;
		cout << "Run with -h to print out the command line options\n" << endl;
		return 0;
	}
	if(b_parameters.errorRate == 0.00 && b_parameters.skipEstimate == true)
	{
		cout << "BELLA execution terminated." 	<< endl;
		cout << " The user should either:" 		<< endl;
		cout << " * -e = suggest an error rate"	<< endl;
		cout << " * -q = confirm that the data has quality values and we can estimate the error rate from the data set" << endl;
		cout << " * -u = confirm that we can use a default error rate (0.15)\n" << endl;
		return 0;
	}
#endif

	free(optList);
	free(thisOpt);
	//
	// Declarations 
	//
	vector<filedata> allfiles = GetFiles(all_inputs_fofn);
	int lower, upper; // reliable range lower and upper bound
	double ratioPhi;
	Kmer::set_k(b_parameters.kmerSize);
	size_t upperlimit = 10000000; // in bytes
	Kmers kmervect;
	vector<string> seqs;
	vector<string> quals;
	vector<string> nametags;
	readVector_ reads;
	Kmers kmersfromreads;
	vector<tuple<size_t,size_t,size_t>> occurrences;
	vector<tuple<size_t,size_t,size_t>> transtuples;
	// 
	// File and setting used
	//
#ifdef PRINT
#ifdef JELLYFISH
	cout << "kmerFile: "				<< kmer_file						<< std::endl;
#endif
	std::cout << std::fixed;
	std::cout << std::setprecision(2);
	std::cout << "outputFile:	"		<< out_file							<< std::endl;
	std::cout << "inputCoverage:	"	<< coverage							<< std::endl;
	std::cout << "kmerSize:	"			<< b_parameters.kmerSize			<< std::endl;
	std::cout << "kmerRift:	"			<< b_parameters.kmerRift			<< std::endl;
	std::cout << "minOverlap:	"		<< b_parameters.minOverlap			<< std::endl;
	std::cout << "minNumKmers:	"		<< b_parameters.minSurvivedKmers	<< std::endl;
	std::cout << "maxOverhang:	"		<< b_parameters.maxOverhang			<< std::endl;
	std::cout << "maxJump:	"			<< b_parameters.maxJump				<< std::endl;
	std::cout << "maxDivergence:	"	<< b_parameters.maxDivergence		<< std::endl;
	std::cout << "outputPaf:	"		<< b_parameters.outputPaf			<< std::endl;
	std::cout << "binSize:	"			<< b_parameters.binSize				<< std::endl;
	std::cout << "deltaChernoff:	"	<< b_parameters.deltaChernoff		<< std::endl;
	std::cout << "runAlignment:	"		<< !b_parameters.skipAlignment		<< std::endl;
	std::cout << "seqAn xDrop:	"		<< b_parameters.xDrop				<< std::endl;
#endif

	//
	// Kmer file parsing, error estimation, reliable bounds computation, and k-mer dictionary creation
	//

	dictionary_t countsreliable;

#ifdef JELLYFISH
	// Reliable bounds computation for Jellyfish using default error rate
	double all = omp_get_wtime();
	lower = computeLower(coverage, b_parameters.errorRate, b_parameters.kmerSize);
	upper = computeUpper(coverage, b_parameters.errorRate, b_parameters.kmerSize);

	std::cout << "errorRate:	"				<< b_parameters.errorRate	<< std::endl;
	std::cout << "kmerFrequencyLowerBound:	"	<< lower					<< std::endl;
	std::cout << "kmerFrequencyUpperBound:	"	<< upper					<< std::endl;

	JellyFishCount(kmer_file, countsreliable, lower, upper);
#else
if(b_parameters.useGerbil)
{
  // Reliable range computation within denovo counting
  cout << "\nRunning with up to " << MAXTHREADS << " threads" << endl;
  double all = omp_get_wtime();
  GerbilDeNovoCount("tempDir", all_inputs_fofn, countsreliable, lower, upper, coverage, upperlimit, b_parameters);
}
else
{ 
	// Reliable range computation within denovo counting
	std::cout << "numThreads:	"				<< MAXTHREADS	<< "\n"		<< std::endl;
	double all = omp_get_wtime();
	DeNovoCount(allfiles, countsreliable, lower, upper, coverage, upperlimit, b_parameters);
}
#ifdef PRINT
	std::cout << "errorRate:	"				<< b_parameters.errorRate	<< std::endl;
	std::cout << "kmerFrequencyLowerBound:	"	<< lower					<< std::endl;
	std::cout << "kmerFrequencyUpperBound:	"	<< upper					<< std::endl;
if(b_parameters.adapThr)
{
	ratioPhi = adaptiveSlope(b_parameters.errorRate);
	std::cout << "adaptiveThreshold constant:	"	<< ratioPhi * (1-b_parameters.deltaChernoff)	<< "\n" << std::endl;
}
else
	std::cout << "userDefinedThreshold:	"	<< b_parameters.defaultThr	<< "\n" << std::endl;
#endif // PRINT
#endif // DENOVO COUNTING

	//
	// Fastq(s) parsing
	//

	double parsefastq = omp_get_wtime();
	size_t read_id = 0; // read_id needs to be global (not just per file)

	vector < vector<tuple<int,int,int>> > alloccurrences(MAXTHREADS);
	vector < vector<tuple<int,int,int>> > alltranstuples(MAXTHREADS);
	vector < readVector_ > allreads(MAXTHREADS);

	for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++)
	{

		ParallelFASTQ *pfq = new ParallelFASTQ();
		pfq->open(itr->filename, false, itr->filesize);

		size_t fillstatus = 1;
		while(fillstatus)
		{ 
			fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
			size_t nreads = seqs.size();

			#pragma omp parallel for
			for(int i=0; i<nreads; i++) 
			{
				// remember that the last valid position is length()-1
				int len = seqs[i].length();

				readType_ temp;
				nametags[i].erase(nametags[i].begin());     // removing "@"
				temp.nametag = nametags[i];
				temp.seq = seqs[i];     // save reads for seeded alignment
				temp.readid = read_id+i;

				allreads[MYTHREAD].push_back(temp);

				for(int j = 0; j <= len - b_parameters.kmerSize; j++)  
				{
					std::string kmerstrfromfastq = seqs[i].substr(j, b_parameters.kmerSize);
					Kmer mykmer(kmerstrfromfastq.c_str(), kmerstrfromfastq.length());
					// remember to use only ::rep() when building kmerdict as well
					Kmer lexsmall = mykmer.rep();

					int idx; // kmer_id
					auto found = countsreliable.find(lexsmall,idx);
					if(found)
					{
						alloccurrences[MYTHREAD].emplace_back(std::make_tuple(read_id+i,idx,j)); // vector<tuple<read_id,kmer_id,kmerpos>>
						alltranstuples[MYTHREAD].emplace_back(std::make_tuple(idx,read_id+i,j)); // transtuples.push_back(col_id,row_id,kmerpos)
					}
				}
			} // for(int i=0; i<nreads; i++)
			//cout << "total number of reads processed so far is " << read_id << endl;
			read_id += nreads;
		} //while(fillstatus) 
		delete pfq;
	} // for all files

	size_t readcount = 0;
	size_t tuplecount = 0;
	for(int t=0; t<MAXTHREADS; ++t)
	{
		readcount += allreads[t].size();
		tuplecount += alloccurrences[t].size();
	}
	reads.resize(readcount);
	occurrences.resize(tuplecount);
	transtuples.resize(tuplecount);

	size_t readssofar = 0;
	size_t tuplesofar = 0;
	for(int t=0; t<MAXTHREADS; ++t)
	{
		copy(allreads[t].begin(), allreads[t].end(), reads.begin()+readssofar);
		readssofar += allreads[t].size();

		copy(alloccurrences[t].begin(), alloccurrences[t].end(), occurrences.begin() + tuplesofar);
		copy(alltranstuples[t].begin(), alltranstuples[t].end(), transtuples.begin() + tuplesofar);
		tuplesofar += alloccurrences[t].size();
	}

	std::sort(reads.begin(), reads.end());   // bool operator in global.h: sort by readid
	std::vector<string>().swap(seqs);        // free memory of seqs  
	std::vector<string>().swap(quals);       // free memory of quals

#ifdef PRINT
	std::cout << "Fastq parsing took:	"		<< omp_get_wtime()	-	parsefastq << "s" << std::endl;
	std::cout << "Total number of reads:	"	<< read_id	<< "\n" << std::endl;
#endif

	//
	// Sparse matrices construction
	//

	double matcreat = omp_get_wtime();

	size_t nkmer = countsreliable.size();
	CSC<size_t, size_t> spmat(occurrences, read_id, nkmer, 
							[] (size_t& p1, size_t& p2) 
							{
								return p1;
							});
	// remove memory of transtuples
	std::vector<tuple<size_t,size_t,size_t>>().swap(occurrences);

	CSC<size_t, size_t> transpmat(transtuples, nkmer, read_id, 
							[] (size_t& p1, size_t& p2) 
							{
								return p1;
							});
	// remove memory of transtuples
	std::vector<tuple<size_t, size_t, size_t>>().swap(transtuples);

#ifdef PRINT
	std::cout << "Sparse matrix construction took:	" << omp_get_wtime()-matcreat << "s\n" << std::endl;
#endif

	//
	// Overlap detection (sparse matrix multiplication) and seed-and-extend alignment
	//
	spmatPtr_ getvaluetype(make_shared<spmatType_>());
	HashSpGEMM(spmat, transpmat, 
		// n-th k-mer positions on read i and on read j
		// AB: not sure if these id1 and id2 are captured correctly, honestly
		[&b_parameters, &reads] (const int& begpH, const int& begpV, const int& id1, const int& id2)
		{
			spmatPtr_ value(make_shared<spmatType_>());

			std::string& read1 = reads[id1].seq;
			std::string& read2 = reads[id2].seq;

			// GG: function in chain.h
			multiop(value, read1, read2, begpH, begpV, b_parameters.kmerSize);
			return value;
		},
		[&b_parameters, &reads] (spmatPtr_& m1, spmatPtr_& m2, const int& id1, const int& id2)
		{
			// GG: after testing correctness, these variables can be removed
			std::string& readname1 = reads[id1].nametag;
			std::string& readname2 = reads[id2].nametag;

			// GG: function in chain.h
			chainop(m1, m2, b_parameters, readname1, readname2);
			return m1;
		},
		reads, getvaluetype, b_parameters.xDrop, out_file, b_parameters, ratioPhi); 

	std::cout << "Total running time:	" << omp_get_wtime()-all << "s\n" << std::endl;
	return 0;
} 
