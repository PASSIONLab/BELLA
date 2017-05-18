#include "CSC.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <ios>
#include <cstdlib> 	
#include <cmath>
#include <algorithm>
#include <vector>

#define KMER_LENGTH 17

using namespace std;

void computeDivergence(vector<double> & dfst, vector<double> & dsnd, map<double, size_t> & divmap) {

	double div, dmin, dmax;
	map<double, size_t>::iterator it;
	double p_ins = 1.090;
	double p_del = 0.955;
	double Lmin, Lmax;

	for (int i = 0; i < dfst.size(); ++i) 
	{
		for(int j = 0; j < dsnd.size(); ++j) 
		{	
			if(dfst[i] <= dsnd[j]) { 
				dmin = dfst[i]; 
				dmax = dsnd[j];
			} else {
				dmax = dfst[i];
				dmin = dsnd[j];
			}

			Lmax = dmin/p_del; // maximum length of the shortest delta given the probability of deletion
			Lmin = dmax/p_ins; // minimum length of the longest delta given the probability of insertion
	
			if(Lmax > Lmin) 
			{
				div = dmax/dmin;
			
				it = divmap.find(div);
				if(it == divmap.end()) {
					divmap.insert(make_pair(div, 1));
				} else it->second++;
			}
		}
	}
}

void cumulativePlot(map<double, size_t> & in) {

	size_t prev = 0;
	ofstream cum2 ("temp.csv");

	for(auto it = in.begin(); it != in.end(); ++it) {

		it->second = prev + it->second;
		prev = it->second;
	}

    if(cum2.is_open()) {
        for(auto it = in.begin(); it != in.end(); ++it) {
 
            cum2 << it->first << " " << it->second << endl; 
        }
    }
    cum2.close();
}

void filterCSV(ifstream & csvin, ofstream & csvout) {

	string div, val, ndiv, nval;
	streampos oldpos;

	if(csvin.is_open()) {
		while(!csvin.eof()) {

			csvin >> div >> val;
			oldpos = csvin.tellg();

			csvin >> ndiv >> nval;

			if(ndiv.compare(div) == 0) {
				continue;
			} else {
				csvout << div << "," << val << endl;
			}
			csvin.seekg(oldpos);
		}
	}
	csvin.close();
	csvout.close();
}


























