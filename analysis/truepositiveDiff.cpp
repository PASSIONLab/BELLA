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
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h> 
#include <memory>
#include <sstream>

using namespace std;

int main (int argc, char* argv[]) 
{
	ifstream truepositive(argv[1]);
	ifstream truepositiveOverT(argv[2]);
	map<pair<string,string>, bool> checkBella;
	map<pair<string,string>, bool>::iterator it;

	int count = 0;
	if(truepositiveOverT.is_open())
	{
		string line;
        while(getline(truepositiveOverT, line))
        {        
        	count++;        
            stringstream lineStream(line);
            string colName, rowName;

            getline(lineStream, colName, '\t');
            getline(lineStream, rowName, '\t');

            it = checkBella.find(make_pair(colName, rowName));
            if(it == checkBella.end())    
                checkBella.insert(make_pair(make_pair(colName, rowName), true));
        }
	}
	truepositiveOverT.close();
	cout << count << endl;
	count = 0;
	ofstream truepositiveDiff;
    truepositiveDiff.open("diff-truepositives-bella.out", std::ofstream::out | std::ofstream::app);

	if(truepositive.is_open())
	{
		string line;
        while(getline(truepositive, line))
        {           
        count++;     
            stringstream lineStream(line);
            string colName, rowName, nkmer, score;
            string colStart, colEnd, colLen, rowStart, rowEnd, rowLen;

            getline(lineStream, colName, '\t');
            getline(lineStream, rowName, '\t');
            getline(lineStream, nkmer, '\t');
            getline(lineStream, score, '\t');
            getline(lineStream, colStart, '\t');
            getline(lineStream, colEnd, '\t');
            getline(lineStream, colLen, '\t');
            getline(lineStream, rowStart, '\t');
            getline(lineStream, rowEnd, '\t');
            getline(lineStream, rowLen, '\t');

            it = checkBella.find(make_pair(colName, rowName));
            if(it == checkBella.end()) 
            	truepositiveDiff << colName << "\t" << rowName << "\t" << nkmer << "\t" << score << "\t" << colStart  << "\t" << colEnd  << "\t" << colEnd  << "\t" << 
                            rowStart  << "\t" << rowEnd  << "\t" << rowLen << endl;  
        }
	}
	truepositiveDiff.close();
	truepositive.close();
	cout << count << endl;

	return 0;
}

