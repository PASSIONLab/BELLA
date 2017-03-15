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
#include <fstream>
#include <istream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <array>
#include <tuple>
#include <cstring>
#include <string.h>
#include <ios>

using namespace std;

typedef std::tuple<int, string, int> searchtuple;

int main (int argc, char* argv[]) {

	ifstream filein (argv[1]);
	ifstream fastafile (argv[2]);
	ofstream log;
	int length;
	int elem;
	int prec = 0, i, j, w, l;
	int num = 0;
	int h = 0;
	std::size_t found;
	std::string kmer;		
	std::string line;
	std::string previousline;
	std::string substr;
	std::vector<pair <int, string>> data;
	int matches = 0;
	std::vector<pair <int, int>> group;
	std::vector<searchtuple> search;
	
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
	cout << "Tuple created." << endl;

	/* Sorting tuple w/ increasing occurrences */

	std::sort(data.begin(), data.end());
	cout << "Tuple sorted." << endl;
	
	/* Creating groups based on occurrence, tuple <occurrence, last position> */
	
	for(w=0; w<data.size()-1; w++) {
		if(data.at(w).first != data.at(w+1).first) {	
			group.push_back(pair <int, int> (data.at(w).first, w+1));	
		}
	}

	cout << "Groups created." << endl;

	/*cout << "Groups created, coordinates <occurrence, last position (start from 1)>:" << endl;
	for(m=0; m<group.size(); m++) {
			
		cout << group.at(m).first << " " 
		<< group.at(m).second << endl;
			
	}*/

	/* Search kmers in fasta file, overlap allowed */
	
	if(fastafile.is_open()) {


		int prec = 0, found, h, num;
		string line;
		for(i=0; i<group.size(); i++) {
			for(j=prec; j<group.at(i).second; j++) {
				getline(fastafile, line);
				fastafile.ignore(line.size(), '\n'); /* I have to go directly to the second line of the fasta file */	
				string prevline("");
				while(getline(fastafile, line)) {
					cout << " lorenzo e' bello " << endl;
	
					string searchon = prevline + line;
					h = 0;
					do {   
						
						found = searchon.find(data.at(j).second, h);
						if(found != std::string::npos) {
							num++;
							h = found+1;
						}
					} while (found != string::npos);

					prevline = line.substr(line.size()-length + 1, length-1);

				}		
				search.push_back(searchtuple (data.at(j).first, data.at(j).second, num));
				fastafile.clear();					  /* I have to clear the eof flag */
				fastafile.seekg(0, fastafile.beg);
				num = 0;
			}
			prec = group.at(i).second;
								
		}
	} else cout << "Unable to open the fasta file.\n\n";

	fastafile.close();

	log.open("log");

	for(int l=0; l<search.size(); l++) {
			
			log << std::get<0>(search.at(l)) << " "
			    << std::get<1>(search.at(l)) << " "
			    << std::get<2>(search.at(l)) << endl;
			
	}

	log.close();
	cout << "\n";
	
return 0;

}
