/*
 * Copyright 2017 Giulia Guidi
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
#include <unordered_map>
#include <sys/stat.h>

#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include "kmercode/Buffer.h"
#include "kmercode/common.h"
#include "kmercode/fq_reader.h"
#include "kmercode/ParallelFASTQ.h"

#define LSIZE 16000
#define KMER_LENGTH 15

using namespace std;

struct filedata {

    char filename[MAX_FILE_PATH];
    size_t filesize;
};

std::vector<filedata>  GetFiles(char *filename)
{
    int64_t totalsize = 0;
    int numfiles = 0;
    vector< filedata > filesview;
    
    filedata fdata;
    ifstream allfiles(filename);
    if(!allfiles.is_open()) {
        cerr << "Could not open " << filename << endl;
        exit(1);
    }
    allfiles.getline(fdata.filename,MAX_FILE_PATH);
    while(!allfiles.eof())
    {
        struct stat st;
        stat(fdata.filename, &st);
        fdata.filesize = st.st_size;
        
        filesview.push_back(fdata);
        cout << filesview.back().filename << " : " << filesview.back().filesize / (1024*1024) << " MB" << endl;
        allfiles.getline(fdata.filename,MAX_FILE_PATH);
        totalsize += fdata.filesize;
        numfiles++;
    }
    return filesview;
}

struct states {

	int a; // number of kmer with occurence 0 
	int b; // number of kmer with occurence 1 
	int c; // number of kmer with occurence bigger than 1 
	
	states(): a(0), b(0), c(0) {
	}
};

struct node {

	node *A;
	node *C;
	node *G;
	node *T;
	unsigned int count;
	unsigned int groupID; // groupID: occurrence after reads generation (pbsim) and kmer generation (jellyfish) 

	node(): A(NULL), C(NULL), G(NULL), T(NULL), count(0), groupID(0) {
	}
};

typedef std::tuple<int, string, int> finalDataset;	 // <groupID, kmer, #matches> 
typedef std::unordered_map<Kmer, int> dictionary;	 //  <k-mer, #kmers> 
typedef std::vector<Kmer> Kmers;

// Function to add a kmer to build the trie 
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

// Function to build the trie 
void buildTrie(node *trieTree, vector<pair <int, string>> dataset, int length) {

	node *trie = trieTree;
	size_t i;

	for(i=0; i<dataset.size(); i++) {
		
		addKmer(trie, dataset.at(i).second.c_str(), length, dataset.at(i).first);

	}
	
	std::cout << "Trie built" << endl;
}

// Function to search matches 
void matchesSearching(FILE *fastafile, node *root, int kmerLength) {
	
	// long lsize = number of char to be read 
	char buffer[LSIZE];
	size_t length = 0; 
	
	while(!feof(fastafile)) {
		size_t bytes = fread(buffer+length, 1, LSIZE-length, fastafile); 
		length += bytes; // length of the matching, used when the end of the buffer occurs during a possible matching 
		size_t offset;

		// offset = position in the buffer
		for(offset = 0; offset < bytes; offset++) {

			if(buffer[offset] == '\n' && buffer[offset] != 'A' && buffer[offset] != 'C' && buffer[offset] != 'G' && buffer[offset] != 'T')
				continue; 
			
			// root iniziatization 
			node *nodeTrie = root;
			int pos;
			int level;
			
			// level = depth position of correct value, pos = depth position, it takes care of non-valid value, e.g. '\n'
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
				// MATCH 
				nodeTrie->count++;
			}
			
			if(offset+pos >= length) {
				// END OF THE BUFFER BEFORE THE EVALUATION OF A MATCH 
				break;
			}
		}
		
		memcpy(buffer, buffer+offset, sizeof(char)*(length-offset)); // memcpy(destination, origin, #bytes) 
		length = length-offset;
	}
	std::cout << "Search finished" << endl;
}

// Recursive function to save <kmer groupID, kmer string, kmer count> 
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
// Function to sort the dataset based on groupID --> NOT NECESSARY
void sortDataset(vector<finalDataset> &kmerData) {
	
	std::sort(kmerData.begin(), kmerData.end());
	std::cout << "Final dataset sorted (based on groupID)" << endl;

}

// Function to create the dictionary 
void dictionaryCreation(dictionary &kmerDictionary, std::vector<pair <int, Kmer>> &data, std::ofstream &fileout) {

	dictionary::iterator it;
	unsigned int count = 0;

	for(int i = 0; i<data.size()-1; i++) {

		if(data.at(i).first < 9 && data.at(i).first > 3)

			it = kmerDictionary.find(data.at(i).second);

			if(it == kmerDictionary.end()) {
				it = kmerDictionary.insert(std::pair<Kmer, int> (data.at(i).second, count)).first;
				count++;
			}

			/* it = kmerDictionary.find(data.at(i).second.twin());

			if(it == kmerDictionary.end()) {
				it = kmerDictionary.insert(std::pair<string, int> (data.at(i).second.twin(), count)).first;
				count++;
			} */

	}
	
	std::cout << "Dictionary created" << endl;

	/* std::cout << "Dictionary created and saved in Dictionary.csv" << endl;
	
	fileout << "k-mer," << "#kmer" << endl;
	for(auto it = kmerDictionary.begin(); it != kmerDictionary.end(); it++) {
			
		fileout << it->first << "," << it->second << endl;
	} */
}

// De-allocation of the tree 
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
	ofstream fileout ("dictionary.txt");
	int length;
	int elem;
	size_t i;
	char *buffer;
	std::string marco;
	std::string line;
	Kmer::set_k(KMER_LENGTH);
	std::vector<pair <int, Kmer>> data;
	std::vector<pair <int, int>> group;
	std::unordered_map <Kmer, int> kmerDictionary;
	std::vector<filedata> allfiles = GetFiles(argv[4]);
    size_t upperlimit = 1000; // 1 thousand reads at a time
    std::vector<string> seqs;
    std::vector<string> quals; // NOT NECESSARY I GUESS
    Kmers kmers;
	// struct node *trieTree;
	// std::vector<finalDataset> kmerData;
	// statesMap statesData;	

	if(argc == 5){
	
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

	// Creating tuple <occurrence, kmer>

	if(filein.is_open()) {
			while(getline(filein, line)) {
		
				if(line.length() == 0)
					break;

				string substring = line.substr(1);
				elem = stoi(substring);
				getline(filein, marco);
				if(elem > 3 && elem < 9) {	
					data.push_back(pair <int, Kmer> (elem, marco); // TO FIX
			
				}										
			}
	} else cout << "Unable to open the input file.\n\n";
	
	filein.close();
	std::cout << "Initial dataset created" << endl;
    
    for(auto itr= allfiles.begin(); itr != allfiles.end(); itr++)
    {
        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);
        
        size_t fillstatus = 1;
        while(fillstatus)
        {
            fillstatus = pfq->fill_block(seqs, quals, upperlimit); // CHECK quals
            size_t nreads = seqs.size();
            
            for(size_t i=0; i<nreads; ++i)
            {
                size_t found;
                found = seqs[i].length();
                
                // Skip this sequence if the length is too short
                if (seqs[i].length() <= KMER_LENGTH) {
                    continue;
                }

                int nkmers = (seqs[i].length()-KMER_LENGTH+1);

                kmers = Kmer::getKmers(seqs[i]); // Calculate all the kmers
                assert(kmers.size() == nkmers);
               
                // Debugging

				std::cout << "Printing the first k-mer and its reverse complement from this read: ";
 		  	    std::cout << kmers[0].toString() << " " << kmers[0].twin().toString() << endl;
            }
            
        }
    }

	// Trie tree allocation 
	// trieTree = (struct node*)calloc(1, sizeof(struct node));

	// Trie tree building 
	// buildTrie(trieTree, data, length);

	// Search matches in trie tree 
	// fastafile = fopen(argv[2], "r");
	// if(fastafile != NULL) {
	//	matchesSearching(fastafile, trieTree, length);
	// } else cout << "Unable to open the fasta file\n";
	
	// fclose(fastafile);	
    // Generation of the final dataset 
	// buffer = (char*)malloc(sizeof(char)*(length+1));
	// buffer[length] = '\0';
	// newDatasetGeneration(trieTree, length, 0, buffer, kmerData);
	// free(buffer); 
	// cout << "Final dataset created" << endl; 	

	// Trie tree de-allocation
	// freeTrie(trieTree);

	// Final dataset sorted 
	//sortDataset(kmerData);

	if(fileout.is_open()) {
		dictionaryCreation(kmerDictionary, data, fileout);
	} else std::cout << "Unable to open the output file\n";
	
	fileout.close();
	cout << "\n";
	
return 0;

}