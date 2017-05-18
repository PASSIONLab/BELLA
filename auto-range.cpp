#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <cstdlib>
#include <algorithm>

#define GENOMESIZE 4641652
#define KMER_LENGTH 17

using namespace std;

int main (int argc, char* argv[]) {

	int rangeStart;
	std::ifstream filein(argv[1]);
	int elem;
	string line;
	map<int,size_t> frequencies;
	map<int,size_t>::iterator it;

	if(filein.is_open()) 
	{
		while(getline(filein, line)) 
		{
			if(line.length() == 0)
				break;

			string substring = line.substr(1);
			elem = stoi(substring);

			getline(filein,line);

			it = frequencies.find(elem);
			if(it == frequencies.end())
			{
				frequencies.insert(make_pair(elem, 1));
			} else it->second++;
		}
	} else std::cout << "unable to open the input file\n";

	filein.close();

	size_t cumSum = 0;
	it = prev(frequencies.end());

	while(cumSum*KMER_LENGTH < GENOMESIZE)
	{	
		cumSum = cumSum + it->second;
		--it;
	}

	rangeStart = it->first;
	cout << "rangeStart = " << rangeStart <<endl;
}