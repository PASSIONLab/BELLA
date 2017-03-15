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
#include <vector>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <array>
#include <tuple>
#include <cstring>
#include <string.h>

using namespace std;
typedef std::vector<string> data;

struct node {

	node *A = NULL;
	node *C = NULL;
	node *G = NULL;
	node *T = NULL;
	unsigned int count = 0;
	boolean fullKmer = false; 
};

void buildTrie(struct node *trieTree, data *dataset, int length) {

	node *trie = trieTree;
	int i;

	for(i=0; i<dataset.size(); i++) {
		
		addKmer(trie, dataset.at(i).c_str(), length);

	}

}

void addKmer(struct node *trieTree, char *kmer, int length) {

	node *traverse = trieTree;
	int i;

	for(i=0; i<length; i++) {
	
		if((kmer[i] == 'A') && (traverse->A == NULL)) {
			traverse->A = (struct node*)calloc(1, sizeof(struct node));
			traverse = traverse->A;
		}
		else if((kmer[i] == 'C') && (traverse->C == NULL)) {
			traverse->C = (struct node*)calloc(1, sizeof(struct node));
			traverse = traverse->C;
		}
		else if((kmer[i] == 'G') && (traverse->G == NULL)) {
			traverse->G = (struct node*)calloc(1, sizeof(struct node));
			traverse = traverse->G;
		}
		else if((kmer[i] == 'T') && (traverse->T == NULL)) {
			traverse->T = (struct node*)calloc(1, sizeof(struct node));
			traverse = traverse->T;
		}
	}
	
	
}

int main () {

	struct node *trieTree = (struct node*)calloc(1, sizeof(struct node));
	char *word;

	cout << "Enter the kmer: ";
	scanf("%s", word);
	buildTree(trieTree, word, 6);

	return 0;
}










