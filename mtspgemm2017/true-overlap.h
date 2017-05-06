#include "CSC.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>

#define KMER_LENGTH 17
#define p_ins 1.09
#define p_del 0.955

using namespace std;

template <typename IT>
bool trueoverlapFilter(vector<IT> & dfst, vector<IT> & dsnd) {

	IT dmin, dmax; // to define the shortest and longest delta
	double Lmin, Lmax; // to compute the shortest length of the longest delta and the longest length of the shortes delta 

	// to discharge those pairs which have an unbalanced presence of k-mers 
	// (i.e. kmers_1 appears 4 times on read i and 2 times on read j)
	if(dfst.size() != dsnd.size()) {
		return false;
	} else {

		// assuming same lengths for the two deltas (given the above if) and assuming positions saved in order
		for (IT i = 0; i < dfst.size(); ++i) 
		{
			if(dfst[i] < dsnd[i]) { 
				dmin = dfst[i]; 
				dmax = dsnd[i];
			} else {
				dmax = dfst[i];
				dmin = dsnd[i];
			}

			Lmax = (double)dmax/(double)(p_del); // maximum length of the shortest delta given the probability of deletion
			Lmin = (double)dmin/(double)(p_ins); // minimum length of the longest delta given the probability of insertion
			cout << Lmax << " Lmax" << endl;
			cout << Lmin << " Lmin" << endl;

			if(Lmax < Lmin) {
				return false;
			}
		}
		return true;
	}
}
