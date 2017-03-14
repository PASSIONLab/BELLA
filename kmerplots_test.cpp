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
	//ofstream fileout ("tuple.txt");
	ifstream fastafile (argv[2]);
	ofstream log;
	int depth;
	int length;
	int elem;
	int prec = 0, k, t, i, j, w, l, num = 0;
	string kmer;		
	string line;
	std::vector<pair <int, string>> data;
	int matches = 0;
	std::vector<pair <int, int>> group;
	std::vector<searchtuple> search;
	
	if(argc == 5){
	
		depth = atoi(argv[3]);
		length = atoi(argv[4]);
		
	}
	else{
		perror("Erroneous input.");
		exit(1);
	}

	cout << "\nThe input file is: " << argv[1] <<endl;
	cout << "The psbsim depth is: " << depth <<endl;
	cout << "The k-mer length is: " << length <<endl;
	cout << "The reference genome is: Escherichia coli, " << argv[2] <<endl;
	cout << "\n";

	/* Creating tuple <occurrence, kmer> */

	if(filein.is_open()) {
		//if(fileout.is_open()) {

			while(getline(filein, line)) {
		
				if(line.length() == 0)
					break;

				string substring = line.substr(1);
				elem = stoi(substring);
				getline(filein, kmer);
				if(elem != 1) {	
					//fileout << elem << "\t" << kmer << endl;
					data.push_back(pair <int, string> (elem, kmer));
			
				}										
			}
		//} else cout << "Unable to open the new output file.\n\n";
	} else cout << "Unable to open the input file.\n\n";
	
	filein.close();
	//fileout.close();

	cout << "Tuple created." << endl;

	/* Sorting tuple w/ increasing occurrences */

	std::sort(data.begin(), data.end());

	/*for(int t = 0; t<data.size(); t++) {
		cout << data.at(t).second << endl;
	}*/

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

	/* Search kmers in fasta file */
	/* To do: fix matches count in reference genome */

	/*if(fastafile.is_open()) {
		while(getline(fastafile, line)) {
			for(k=0; k<=(line.size()-length); k++) {

				std::string substr = line.substr(k,length);
				//const char *substr_cr = substr.c_str();
      				if((strcmp(substr.c_str(), "GCATT") == 0))
					cout << "GCATT FOUND" <<endl;
							
			}
		}		
	} else cout << "Unable to open the fasta file.\n\n";*/

	/* To fix: data.at(j).second is not updated into the while loop 
	   To do: implement k-mer search between two lines */
	
	/* cout << data.size() << endl; --> 770 */
	/*cout << group.size() << endl; --> 46 */

	if(fastafile.is_open()) {
		for(i=0; i<group.size(); i++) {	
			for(j=prec; j<group.at(i).second; j++) {
				while(getline(fastafile, line)) {
					//cout << data.at(j).second.c_str() << endl;
					for(k=0; k<=(line.size()-length); k++) {

						std::string substr = line.substr(k,length);
						//const char *substr_cr = substr.c_str();
						if((strcmp(substr.c_str(), data.at(j).second.c_str())) == 0)
							num++;
					//	cout << num << " " << substr.c_str() << " " << data.at(j).second.c_str() << endl;
							
					}
				}
				prec = group.at(i).second;		
				search.push_back(searchtuple (data.at(j).first, data.at(j).second, num));
				num = 0;
				fastafile.clear();
				fastafile.seekg(0, fastafile.beg);
				/*cout << j << " " << group.at(i).second << endl; /* out of range */
			}					
		}
	} else cout << "Unable to open the fasta file.\n\n";

	fastafile.close();

	/*cout << "Length of search: "  << search.size() << endl;*/

	log.open("log");

	for(l=0; l<search.size(); l++) {
			
			log << std::get<0>(search.at(l)) << " "
			    << std::get<1>(search.at(l)) << " "
			    << std::get<2>(search.at(l)) << endl;
			
	}

	log.close();
	cout << "\n";
	
return 0;

}
