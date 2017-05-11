#include "CSC.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <vector>

#define KMER_LENGTH 17
#define p_ins 1.09
#define p_del 0.955

using namespace std;

bool trueoverlapFilter(vector<double> & dfst, vector<double> & dsnd) {

	double dmin, dmax; // to define the shortest and longest delta
	double Lmin, Lmax; // to compute the shortest length of the longest delta and the longest length of the shortes delta 
	bool same_region = false;

	cout << "dfst.size() = " << dfst.size() << endl;
	cout << "dsnd.size() = " << dsnd.size() << endl;

	for (int i = 0; i < dfst.size(); ++i) 
	{
		for(int j = 0; j < dsnd.size(); ++j)
		{
			if(dfst[i] = dsnd[j]) {
				same_region = true;
			} else {

				if(dfst[i] < dsnd[j]) { 
					dmin = dfst[i]; 
					dmax = dsnd[j];
				} else {
					dmax = dfst[i];
					dmin = dsnd[j];
				}

				Lmax = dmin/p_del; // maximum length of the shortest delta given the probability of deletion
				Lmin = dmax/p_ins; // minimum length of the longest delta given the probability of insertion

				if(Lmax > Lmin) {
					same_region = true;
				}
			}
		}
	}
	return same_region;
}

