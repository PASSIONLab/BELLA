#include "CSC.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>

#define KMER_LENGTH 17

using namespace std;

template <typename IT>
vector<IT> computeDelta(vector<IT> & fst, vector<IT> & snd) {

    vector<IT> delta;
    
    #pragma omp parallel for
    for(IT i = 0; i < fst.size(); ++i) 
    { 
        for(IT j = 0; j < snd.size(); ++j)
        {
            delta.push_back(snd[j]-fst[i]); 
        }
    }

    return delta;
}            
