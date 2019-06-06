/**
* Evaluation code to compare long read-to-read aligner recall/precision
* @author: Giulia Guidi
*/

#ifdef __cplusplus
extern "C" {
#endif

#include "../optlist/optlist.h" // command line parser

#ifdef __cplusplus
}
#endif

#include "IntervalTree.h"
#include "evaluation.h"
#include "def.h"
#include <omp.h>
#include <set>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <sys/types.h> 
#include <sys/stat.h> 
#include <math.h>
#include <limits.h>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h> 
#include <memory>

int main (int argc, char* argv[]) {

	option_t *optList, *thisOpt;
	optList = NULL;
	optList = GetOptList(argc, argv, (char*)"G:B:m:M:zl:i:H:D:L:");

	int  	minOverlap  = 2000;
	bool 	isSimulated = false;
	char*	G = NULL;
	char*	B = NULL; 	// BELLA standard output
	char*	m = NULL; 	// minimap/minimap2
	char*	M = NULL; 	// MECAT
	char*	i = NULL; 	// MECAT's indexes
	char*	H = NULL; 	// MHAP
	char*	D = NULL; 	// DALIGNER
	char*	L = NULL; 	// BLASR

	if(optList == NULL) {
		std::cout << "Program execution terminated: not enough parameters or invalid option." 	<< std::endl;
		std::cout << "Run with -h to print out the command line options.\n" 					<< std::endl;
		return 0; }

	while (optList!=NULL) {
		thisOpt = optList;
		optList = optList->next;
		switch (thisOpt->option) {
			case 'G': {
				if(thisOpt->argument == NULL)
				{
					std::cout << "\nProgram execution terminated: -G requires the ground truth file." 	<< std::endl;
					std::cout << "Run with -h to print out the command line options.\n" 				<< std::endl;
					return 0;
				}
				G = strdup(thisOpt->argument);
				break;
			}
			case 'l': {
				minOverlap = stoi(thisOpt->argument);
				break;
			}
			case 'z': {
				isSimulated = true;
				break;
			}
			case 'B': {
				B = strdup(thisOpt->argument);
				break;
			}
			case 'm': {
				m = strdup(thisOpt->argument);
				break;
			}
			case 'M': {
				M = strdup(thisOpt->argument);
				break;
			}
			case 'i': {
				i = strdup(thisOpt->argument);
				break;
			}
			case 'H': {
				H = strdup(thisOpt->argument);
				break;
			}
			case 'L': {
				L = strdup(thisOpt->argument);
				break;
			}
			case 'D': {
				D = strdup(thisOpt->argument);
				break;
			}
			case 'h': {
				std::cout << "\nUsage:\n" 							<< std::endl;
				std::cout << " -z : Simulated data [false]" 		<< std::endl;
				std::cout << " -l : minOverlap length [2000]" 		<< std::endl;
				std::cout << " -G : Ground truth file (required)" 	<< std::endl;
				std::cout << " -B : BELLA standard output format" 	<< std::endl;
				std::cout << " -m : minimap/minimap2/PAF format" 	<< std::endl;
				std::cout << " -H : MHAP output format" 			<< std::endl;
				std::cout << " -D : DALIGNER output format" 		<< std::endl;
				std::cout << " -L : BLASR output format" 			<< std::endl;
				std::cout << " -M : MECAT output format" 			<< std::endl;
				std::cout << " -i : MECAT's indexes\n" 				<< std::endl;
				FreeOptList(thisOpt);
				return 0;
			}
		}
	}

	if(G == NULL) {
		cout << "\nProgram execution terminated: missing argument." 	<< endl;
		cout << "Run with -h to print out the command line options.\n" 	<< endl;
		return 0;
	}

	if(M != NULL && i == NULL) {
		cout << "\nProgram execution terminated: missing argument." 													<< std::endl;
		cout << "Please add MECAT idx2read file with -i option. Run with -h to print out the command line options.\n" 	<< std::endl;
		return 0;
	}

	free(optList);
	free(thisOpt);

	bool duplicate = false; // some software doesn't output both (A,B) and (B,A)
	std::ifstream data(G);
	std::multiset<entry, classcom> Gset = readTruthOutput(data, minOverlap, isSimulated); // ground truth
	std::multiset<entry, classcom> Sset; // software output

	if(B) {
		std::ifstream reads(B);
		Sset = readBellaOutput(reads);
		duplicate = true;
		std::cout << "Bella" << std::endl;
		evaluate(Sset, Gset, minOverlap, duplicate);
	}

	if(m) {
		std::ifstream reads(m);
		Sset = readMinimapOutput(reads);
		duplicate = true;
		std::cout << "Minimap2" << std::endl;
		evaluate(Sset, Gset, minOverlap, duplicate);
	}

	if(M) {
		std::ifstream reads(M);
		std::ifstream index(i);
		Sset = readMecatOutput(reads, index);
		duplicate = true;
		std::cout << "Mecat" << std::endl;
		evaluate(Sset, Gset, minOverlap, duplicate);
	}

	if(H) {
		std::ifstream reads(H);
		Sset = readMhapOutput(reads);
		duplicate = false;
		std::cout << "Mhap" << std::endl;
		evaluate(Sset, Gset, minOverlap, duplicate);
	}

	if(L) {
		std::ifstream reads(L);
		Sset = readBlasrOutput(reads);
		duplicate = false;
		std::cout << "Blasr" << std::endl;
		evaluate(Sset, Gset, minOverlap, duplicate);
	}

	if(D) {
		std::ifstream reads(D);
		Sset = readDalignerOutput(reads);
		duplicate = false;
		std::cout << "Daligner" << std::endl;
		evaluate(Sset, Gset, minOverlap, duplicate);
	}
}
