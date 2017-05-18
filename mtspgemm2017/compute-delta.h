#include "CSC.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <ctgmath>
#include <algorithm>

#define KMER_LENGTH 17

using namespace std;

template <typename IT>
vector<double> computeDelta(vector<IT> & fst, vector<IT> & snd) {

    vector<double> delta;
    double min, max, diff;

    //#pragma omp parallel for
    for(size_t i = 0; i < fst.size(); ++i) 
    { 
        for(size_t j = 0; j < snd.size(); ++j)
        {
            //cout << "i =" << fst[i] << ", j = " << snd[j] << endl;

            if(fst[i] < snd[j]) {
                min = (double)fst[i];
                max = (double)snd[j];
            } else {
                max = (double)fst[i];
                min = (double)snd[j];
            }

            diff = max-min;
            delta.push_back(diff); 
        }
    }

    return delta;
}            
