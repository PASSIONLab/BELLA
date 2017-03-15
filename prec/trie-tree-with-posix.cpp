/*
 * Copyright 2017 Giulia Guidi, Marco Rabozzi, Alberto Scolari
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

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
#include <stack>
#include <cstring>
#include <string.h>
#include <cassert>
#include <ios>
#include <map>

#define LSIZE 16000
#define GENOMESIZE 4699679

using namespace std;

struct states {

	int a; /* number of kmer with occurence 0 */
	int b; /* number of kmer with occurence 1 */
	int c; /* number of kmer with occurence bigger than 1 */
	
	states(): a(0), b(0), c(0) {
	}
};

struct node {

	node *A;
	node *C;
	node *G;
	node *T;
	unsigned int count;
	unsigned int groupID; /* groupID: occurrence after reads generation (pbsim) and kmer generation (jellyfish) */

	node(): A(NULL), C(NULL), G(NULL), T(NULL), count(0), groupID(0) {
	}
};

typedef std::tuple<int, string, int> finalDataset; /* <groupID, kmer, #matches> */
typedef std::map<int, states> statesMap;           /* <key value = groupID, mapped value = possibleStates> */
//typedef std::vector<int> posinGenome (4699731);  /* vector containing all the genome positions */
typedef std::map<int, std::vector<int>> groupPos;  /* <group, where the matches are in the genome> */

/* Function to add a kmer to build the trie */
void addKmer(node *trieTree, const char *kmer, int length, int groupNum) {

	node *traverse = trieTree;
	int i;

	for(i=0; i<length; i++) {
	
		if(kmer[i] == 'A') {
			if(traverse->A == NULL) {
				traverse->A = (struct node*)calloc(1, sizeof(struct node));
			}
			traverse = traverse->A;
		}	
		else if(kmer[i] == 'C') {
			if(traverse->C == NULL) {
				traverse->C = (struct node*)calloc(1, sizeof(struct node));
			}
			traverse = traverse->C;
		}
		else if(kmer[i] == 'G') {
			if(traverse->G == NULL) {
				traverse->G = (struct node*)calloc(1, sizeof(struct node));
			}
			traverse = traverse->G;
		}
		else if(kmer[i] == 'T') {
			if(traverse->T == NULL) {
				traverse->T = (struct node*)calloc(1, sizeof(struct node));
			}
			traverse = traverse->T;
		}
	}
	traverse->groupID = groupNum;
}

/* Function to build the trie */
void buildTrie(node *trieTree, vector<pair <int, string>> dataset, int length) {

	node *trie = trieTree;
	size_t i;

	for(i=0; i<dataset.size(); i++) {
		
		addKmer(trie, dataset.at(i).second.c_str(), length, dataset.at(i).first);

	}
	
	std::cout << "Trie built" << endl;
}

/* Function to search matches */
void matchesSearching(FILE *fastafile, node *root, int kmerLength, groupPos &dataPos, std::ofstream &filepos) {
	
	/* long lsize = number of char to be read */
	char buffer[LSIZE];
	size_t length = 0; 
	size_t posix = 0; /* to track where in the genome the match occurs */
	groupPos::iterator it;	

	while(!feof(fastafile)) {
		size_t bytes = fread(buffer+length, 1, LSIZE-length, fastafile); 
		length += bytes; /* length of the matching, used when the end of the buffer occurs during a possible matching */
		size_t offset;

		/* offset = position in the buffer */
		for(offset = 0; offset < bytes; offset++) {
		posix++;
		//std::cout << posix << endl;

			if(buffer[offset] == '\n' && buffer[offset] != 'A' && buffer[offset] != 'C' && buffer[offset] != 'G' && buffer[offset] != 'T')
				continue; 
			
			/* root iniziatization */
			node *nodeTrie = root;
			int pos;
			int level;
			
			/* level = depth position of correct value, pos = depth position, it takes care of non-valid value, e.g. '\n' */
			for(level = 0, pos = 0; (level<kmerLength) && (nodeTrie != NULL) && ((offset+pos)<length); pos++) {
			
				if(buffer[offset+pos] == '\n' && buffer[offset+pos] != 'A' && buffer[offset+pos] != 'C' && buffer[offset+pos] != 'G' && buffer[offset+pos] != 'T')
					continue;
			
				level++;
				
				char base = buffer[offset+pos];
				
				if(base == 'A')
					nodeTrie = nodeTrie->A;
				else if(base == 'C')
					nodeTrie = nodeTrie->C;
				else if(base == 'G')
					nodeTrie = nodeTrie->G;
				else if(base == 'T')
					nodeTrie = nodeTrie->T;
				else break;

			}
			
			if(nodeTrie != NULL && level == kmerLength) {	
				/* MATCH */
				nodeTrie->count++;
				it = dataPos.find(nodeTrie->groupID);
			
				if(it == dataPos.end()) {
				it = dataPos.insert(std::pair<int, std::vector<int>> (nodeTrie->groupID, std::vector<int>())).first; /* Chiamo il costruttore dell'oggetto che mi serve */
				}

				it->second.push_back(posix-kmerLength+1);
			}
			
			if(offset+pos >= length) {
				/* END OF THE BUFFER BEFORE THE EVALUATION OF A MATCH */
				break;
			}
		}
		
		memcpy(buffer, buffer+offset, sizeof(char)*(length-offset)); /* memcpy(destination, origin, #bytes) */
		length = length-offset;
	}
	std::cout << "Search finished and matches reported in posinGenome.csv" << endl;

	filepos << "GroupID," << "Position of the match in the genome" << endl;
				filepos << "\n" << endl;
				for(auto it = dataPos.begin(); it != dataPos.end(); it++) {
			
					filepos << it->first << ","; 
					for(size_t p = 0; p<it->second.size(); p++) {
						filepos << it->second[p] << ","; 
					}
						filepos << endl;
				}
}

/* Recursive function to save <kmer groupID, kmer string, kmer count> */
void newDatasetGeneration(node *trieTree, int kmerLength, int level, char *buffer, vector<finalDataset> &kmerData) {

	if(level == kmerLength) {
		kmerData.push_back(finalDataset (trieTree->groupID, string(buffer), trieTree->count));
		return;
	}

	if(trieTree->A != NULL) {	
		buffer[level] = 'A';
		newDatasetGeneration(trieTree->A, kmerLength, level+1, buffer, kmerData);
	} 
	
	if(trieTree->C != NULL) {
		buffer[level] = 'C';
		newDatasetGeneration(trieTree->C, kmerLength, level+1, buffer, kmerData);
	}

	if(trieTree->G != NULL) {
		buffer[level] = 'G';
		newDatasetGeneration(trieTree->G, kmerLength, level+1, buffer, kmerData);
	}
	
	if(trieTree->T != NULL) {
		buffer[level] = 'T';
		newDatasetGeneration(trieTree->T, kmerLength, level+1, buffer, kmerData);
	} 
}
/* Function to sort the dataset based on groupID --> NOT NECESSARY */
void sortDataset(vector<finalDataset> &kmerData) {
	
	std::sort(kmerData.begin(), kmerData.end());
	std::cout << "Final dataset sorted (based on groupID)" << endl;

}

/* Function to divide the dataset into the three different possible states */
void statesDivision(vector<finalDataset> &kmerData, statesMap &statesData, std::ofstream &fileout) {

	statesMap::iterator it;

	for(int i = 0; i<kmerData.size()-1; i++) {

		it = statesData.find(std::get<0>(kmerData.at(i)));

		if(it == statesData.end()) {
			states groupStates;
			it = statesData.insert(std::pair<int, states> (std::get<0>(kmerData.at(i)), groupStates)).first;
		}
	
		switch (std::get<2>(kmerData.at(i))) {
		case 0:
			it->second.a++;
			break;
		case 1:
			it->second.b++;
			break;
		default:
			it->second.c++;
		}	
	}
	std::cout << "States map created and saved in finalDataset.csv" << endl;
	
	fileout << "GroupID," << "StateA," << "StateB," << "StateC" << endl;
	fileout << "\n" << endl;
	for(auto it = statesData.begin(); it != statesData.end(); it++) {
			
		fileout << it->first << "," << it->second.a << "," << it->second.b << "," << it->second.c << endl;
	}
}

/* De-allocation of the tree */
void freeTrie(node *trieTree) {

	if(trieTree->A != NULL) {	
		freeTrie(trieTree->A);
	} 
	
	if(trieTree->C != NULL) {
		freeTrie(trieTree->C);
	}

	if(trieTree->G != NULL) {
		freeTrie(trieTree->G);
	}
	
	if(trieTree->T != NULL) {
		freeTrie(trieTree->T);
	} 

	free(trieTree);
}

int main (int argc, char* argv[]) {

	ifstream filein (argv[1]);
	FILE *fastafile;
	ofstream fileout ("finalDataset.csv");
	ofstream filepos ("posinGenome.csv");
	int length;
	int elem;
	size_t i;
	char *buffer;
	std::string kmer;
	std::string line;
	std::vector<pair <int, string>> data;
	std::vector<pair <int, int>> group;
	struct node *trieTree;
	std::vector<finalDataset> kmerData;
	statesMap statesData;
	groupPos dataPos;

	if(argc == 4){
	
		length = atoi(argv[3]);
		
	}
	else{
		perror("Erroneous input.");
		exit(1);
	}

	cout << "\nThe input file is: " << argv[1] <<endl;
	cout << "The psbsim depth is: 30" << endl;
	cout << "The k-mer length is: " << length <<endl;
	cout << "The reference genome is: Escherichia coli, " << argv[2] <<endl;
	cout << "\n";

	/* Creating tuple <occurrence, kmer> */
	if(filein.is_open()) {
			while(getline(filein, line)) {
		
				if(line.length() == 0)
					break;

				string substring = line.substr(1);
				elem = stoi(substring);
				getline(filein, kmer);
				if(elem != 1) {	
					data.push_back(pair <int, string> (elem, kmer));
			
				}										
			}
	} else cout << "Unable to open the input file.\n\n";
	
	filein.close();
	cout << "Initial dataset created" << endl;

	/* Trie tree allocation */
	trieTree = (struct node*)calloc(1, sizeof(struct node));

	/* Trie tree building */
	buildTrie(trieTree, data, length);

	/* Search matches in trie tree */
	fastafile = fopen(argv[2], "r");
	if((fastafile != NULL) && (filepos.is_open())) {
		matchesSearching(fastafile, trieTree, length, dataPos, filepos);
	} else std::cout << "Unable to open the fasta file or the pos file\n";
	
	fclose(fastafile);	

	/* Generation of the final dataset */
	buffer = (char*)malloc(sizeof(char)*(length+1));
	buffer[length] = '\0';
	newDatasetGeneration(trieTree, length, 0, buffer, kmerData);
	free(buffer); 
	cout << "Final dataset created" << endl;	

	/* Trie tree de-allocation */
	freeTrie(trieTree);

	/* Final dataset sorted */
	//sortDataset(kmerData);

	/* TO DO: division in states */
	if(fileout.is_open()) {
		statesDivision(kmerData, statesData, fileout);
	} else std::cout << "Unable to open the output file\n";
	
	fileout.close();
	cout << "\n";
	
return 0;

}
