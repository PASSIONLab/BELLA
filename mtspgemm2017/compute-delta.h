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
    for(IT i = 0; i < fst.size(); ++i) 
    { 
        for(IT j = 0; j < snd.size(); ++j)
        {
            if(fst[i] < snd[j]) {
                min = fst[i];
                max = snd[j];
            } else {
                max = fst[i];
                min = snd[j];
            }

            diff = max-min;
            delta.push_back(diff); 
        }
    }
    //cout << "delta.size() = " << delta.size() << endl;
    return delta;
}            
