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
#include <map>

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

struct node {

	node *A;
	node *C;
	node *G;
	node *T;
	unsigned int kmerindex;

	node(): A(NULL), C(NULL), G(NULL), T(NULL), kmerindex(0) {
	}
};

typedef std::tuple<int, string, int> finalDataset;	 // <groupID, kmer, #matches> 
typedef std::map<Kmer, int> dictionary;	 			 // <k-mer && reverse-complement, #kmers> 
typedef std::vector<Kmer> Kmers;
typedef std::vector<pair<int,unsigned int>> occurrences;
typedef std::map<string,int> readtoindex_map;

// Function to add a kmer to build the trie 
void addKmer(node *trieTree, const char *kmer, int count) {

	node *traverse = trieTree;
	int i;

	for(i=0; i<KMER_LENGTH; i++) {
	
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
	traverse->kmerindex = count;
}

// Function to build the trie 
void buildTrie(node *trieTree, dictionary &src) {

	node *trie = trieTree;
	dictionary::iterator it;

	for(it=src.begin(); it!=src.end(); it++) {
		
		addKmer(trie, it->first.toString().c_str(), it->second);

	}
}

//Function to search matches between reads and dictionary 
occurrences matchingFoo(std::vector<string> &reads, node *root, readtoindex_map &readtoindex) {

	occurrences indices;
	readtoindex_map::iterator it = readtoindex.begin();

	for (int i=0; i<reads.size(); ++i) {

        Kmers kmersfromreads = Kmer::getKmers(reads[i]);        // calculate all the kmers
        Kmers twin_kmersfromreads;   
        for(int k = 0; k<kmersfromreads.size(); k++) {	 	   
    	    twin_kmersfromreads[k] = kmersfromreads[k].twin(); // calculate all the kmers reverse complement, not sure if I need reverse complement of the reads
        }      
        kmersfromreads.insert(kmersfromreads.end(), twin_kmersfromreads.begin(), twin_kmersfromreads.end()); // append the 2nd vct to the 1st
        int nkmers = (reads[i].length()-KMER_LENGTH+1);
        assert(kmersfromreads.size() == nkmers);

        for(int j=0; j<nkmers; j++) {

        	node *nodeTrie = root;
			int level;
			
			// level = depth position of correct value
			for(level = 0; (level<KMER_LENGTH) && (nodeTrie != NULL); level++) {
				
				const char *base = kmersfromreads[j].toString().c_str();
				
				if(base[level] == 'A')
					nodeTrie = nodeTrie->A;
				else if(base[level] == 'C')
					nodeTrie = nodeTrie->C;
				else if(base[level] == 'G')
					nodeTrie = nodeTrie->G;
				else if(base[level] == 'T')
					nodeTrie = nodeTrie->T;
				else break;

			}
			// MATCH 
			if(nodeTrie != NULL && level == KMER_LENGTH)
				indices.push_back(pair<int,unsigned int> (i, nodeTrie->kmerindex));	
				// here I'm saving the index of the read and the index ok the kmer in dictionary
				// not sure if it's correct --> alternative case: indices.push_back(indexof(i), indexof(kmersfromstring))
		}
		it = readtoindex.find(reads[i]);
		if(it == readtoindex.end())
			it = readtoindex.insert(std::pair<string, int> (reads[i], i)).first; // here I'm saving read-to-index map
	}
	return indices;
}

// Function to create the dictionary 
void dictionaryCreation(dictionary &kmerdict, Kmers &kmervect) {

	dictionary::iterator it;
	unsigned int count = 0;

	for(int i = 0; i<kmervect.size()-1; i++) { // TO FIX: EXIT THE LOOP BUT TAKE A LOT OF TIME
		it = kmerdict.find(kmervect[i]);

		if(it == kmerdict.end()) {
			it = kmerdict.insert(std::pair<Kmer, int> (kmervect[i], count)).first;
			count++;
		}
	}
	std::cout << "kmervect in the dict" << endl;

	for(int i = 0; i<kmervect.size()-1; i++) {

		it = kmerdict.find(kmervect[i].twin());

		if(it == kmerdict.end()) {
			it = kmerdict.insert(std::pair<Kmer, int> (kmervect[i].twin(), count)).first;
			count++;
		}
	}
	std::cout << "kmervect.twin() in the dict" << endl;
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
	ofstream fileout ("out.csv");
	int elem;
	size_t i;
	char *buffer;
	std::string kmerstr;
	std::string line;
	Kmer::set_k(KMER_LENGTH);
	dictionary kmerdict;
	std::vector<filedata> allfiles = GetFiles(argv[3]);
    size_t upperlimit = 1000; // 1 thousand reads at a time
    Kmer kmerfromstr;
    Kmers kmervect;
    std::vector<string> seqs;
    std::vector<string> quals;
    Kmers kmersfromreads;
	struct node *trieTree; 

	cout << "\ninput file : " << argv[1] <<endl;
	cout << "psbsim depth : 30" << endl;
	cout << "k-mer length : " << KMER_LENGTH <<endl;
	cout << "reference genome : escherichia coli, " << argv[2] <<endl;
	// filtering on kmers --> kmers which occur between 4 and 8 times 
	if(filein.is_open()) {
			while(getline(filein, line)) {
				if(line.length() == 0)
					break;

				string substring = line.substr(1);
				elem = stoi(substring);
				getline(filein, kmerstr); 
				if(elem>3 && elem<9) {	
					//strtoChar = kmerstr.c_str();
					kmerfromstr.set_kmer(kmerstr.c_str());
					kmervect.push_back(kmerfromstr);
				}									
			}
	} else cout << "Unable to open the input file.\n\n";
	filein.close();

	std::cout << "filtered dataset parsed of size: " << kmervect.size() << " elem" << endl;
	dictionaryCreation(kmerdict, kmervect);

    for(auto itr= allfiles.begin(); itr != allfiles.end(); itr++) {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);
        
        size_t fillstatus = 1;
        while(fillstatus) {

            fillstatus = pfq->fill_block(seqs, quals, upperlimit);
            size_t nreads = seqs.size();
            
            for(size_t i=0; i<nreads; ++i) {
                size_t found;
                found = seqs[i].length();
                
                // skip this sequence if the length is too short
                if (seqs[i].length() <= KMER_LENGTH)
                    continue;

                // debugging
    			// std::cout << seqs.size() << endl; // always 1
            }
            // std::cout << seqs.size() << endl; always 1 until the last iteration which is equal to 0 
        } 
    }
    std::cout << "fastq file parsed" << endl;
    // debugging
    if(seqs.empty()) 
    	std::cout << "seqs vector is empty" << endl; // it's an empty vector

	//trie tree allocation 
	trieTree = (struct node*)calloc(1, sizeof(struct node));

	//trie tree building 
	buildTrie(trieTree, kmerdict);
	std::cout << "kmers trie built" << endl;

	//search matching
	readtoindex_map readtoindex;
	occurrences matches = matchingFoo(seqs, trieTree, readtoindex);
	std::cout << "search ended : occurrences vector and read-to-index map created" << endl;

	/*if(fileout.is_open()) {
		for(int s = 0; s<matches.size(); s++) {
			fileout << matches[s].first << "," << matches[s].second << endl; // doesn't work --> seqs is an empty vector
		}
	}*/

	freeTrie(trieTree);

return 0;
}