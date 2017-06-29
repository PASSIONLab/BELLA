#include "CSC.h"
#include "utility.h"
#include "IntervalTree.h"
#include "BitMap.h"
#include "global.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#define KMER_LENGTH 17
#define S 1 // minimum number of shard k-mers
using namespace std;

// compute the number of true overlapping reads
double TrueOverlapping(ifstream & ifs) {

    vector<Interval<int>> intervals;
    vector<Interval<int>> queries;
    vector<Interval<int>>::iterator q;
    double trueoverlaps;

    if(ifs.is_open()) {

        int id;
        int st;
        int en;
        
        while(ifs >> id >> st >> en) {
            
            intervals.push_back(Interval<int>(st, en, id));
            queries.push_back(Interval<int>(st, en, id));
        }

    } else std::cout << "error opening the ifs" << endl;

    IntervalTree<int> tree;
    vector<size_t> treecounts;

    tree = IntervalTree<int>(intervals);

    for (q = queries.begin(); q != queries.end(); ++q) 
    {
        vector<Interval<int>> results;
        tree.findOverlapping(q->start, q->stop, results);
        treecounts.push_back(results.size());
    }

    for(size_t t = 0; t < treecounts.size(); ++t) {
        trueoverlaps = trueoverlaps + (double)treecounts[t];
    }

    return trueoverlaps;

}

// compute the distance between two k-mers on the same read
template <typename IT>
double ComputeDistance(IT & fst, IT & snd) {

    double delta;
    double min, max, diff;

    //#pragma omp parallel for
    //cout << "i =" << fst[i] << ", j = " << snd[j] << endl;

    if(fst < snd) {
        min = (double)fst;
        max = (double)snd;
    } else {
        max = (double)fst;
        min = (double)snd;
    }

    diff = max-min;
    delta = diff; 

    return delta;
}    

// find when a k-mers pair represents an evidence of potential overlap or not
bool PotentialOverlap(double & dfst, double & dsnd) {

    double dmin, dmax; // to define the shortest and longest delta
    double Lmin, Lmax; // to compute the shortest length of the longest delta and the longest length of the shortes delta 
    bool region = false;
    double p_ins = 1.090;
    double p_del = 0.955;
    
    if(dfst <= dsnd) { 
        dmin = dfst; 
        dmax = dsnd;
    } else {
        dmax = dfst;
        dmin = dsnd;
    }
    
    Lmax = dmin/p_del; // maximum length of the shortest delta given the probability of deletion
    Lmin = dmax/p_ins; // minimum length of the longest delta given the probability of insertion
    
    if(Lmax > Lmin) {
        region = true;
    } else region = false;

    return region;
}

// compute the overlap length between potential overlapping reads pairs
template <typename IT>
size_t ComputeLength(map<size_t, pair<size_t, size_t>> & ifsmap, IT & col, IT & row) {

    size_t alignment_length = 0;
    std::map<size_t, pair<size_t, size_t>>::iterator jit;
    std::map<size_t, pair<size_t, size_t>>::iterator iit;

    jit = ifsmap.find(col); // col index
    iit = ifsmap.find(row); // row index
    
    if(iit->second.first < jit->second.first) {
        if(iit->second.second > jit->second.first) {
            alignment_length = min((iit->second.second - jit->second.first), (jit->second.second - jit->second.first));
        }
    }
    else if (iit->second.first > jit->second.first) {
        if(jit->second.second > iit->second.first) {
            alignment_length = min((jit->second.second - iit->second.first), (iit->second.second - iit->second.first));
        }
    } else { 
        alignment_length = min((jit->second.second - iit->second.first), (iit->second.second - iit->second.first)); 
    } 

    return alignment_length;
}

// Estimation of overlapping region length
template <typename IT>
double ExpectedKmers(IT & i, IT & j, IT & len_i, IT & len_j) { // reads length has to be included and position of the same k-mer on different reads

    double left;    // to define the left margin 
    double right;   // to define the right margin
    double estime, p;  // expected overlap region length
    size_t temp_i, temp_j;
    double er = 0.15; // error rate
    
    if(i <= j)
        left = (double)i;
    else left = (double)j;

    temp_i = len_i - i;
    temp_j = len_j - j; 

    if(temp_i <= temp_j)
        right = (double)temp_i;
    else right = (double)temp_j;

    estime = left + right; // estimated overlap
    p = 1-pow((1-pow((1-er), (2*KMER_LENGTH))), estime); // expected number of k-mers

    return p;
}

template <typename IT>
double MaxGap(IT & i, IT & j, IT & len_i, IT & len_j) { // reads length has to be included and position of the same k-mer on different reads

    double left;    // to define the left margin 
    double right;   // to define the right margin
    double estime;  // expected overlap region length
    double gap;
    size_t temp_i, temp_j;
    double variance_single_base = 0.12;
    
    if(i <= j)
        left = (double)i;
    else left = (double)j;

    temp_i = len_i - i;
    temp_j = len_j - j; 

    if(temp_i <= temp_j)
        right = (double)temp_i;
    else right = (double)temp_j;

    estime = left + right; // estimated overlap
    //gap = estime*variance_single_base;
    
    return estime;
}

template <class IT, class NT>
void GetStatistics(const CSC<IT,NT> & A) 
{
    std::ifstream ifs("test_01.axt"); // it would be better to make this reading more general
    std::map<size_t, pair<size_t, size_t>> ifsmap;
    double di; // k-mers distances on read i
    double dj; // k-mers distances on read j
    double TPandFP = 0, TP = 0;
    size_t track = 0, var, alignment_length;
    bool same;

    // creation reads map from axt file to be used to find potential overlapping reads pairs 
    if(ifs.is_open()) {
        int id;
        int st;
        int en;
        std::pair<size_t, size_t> values;
        while(ifs >> id >> st >> en) {
            size_t key = {id};
            size_t start = {st};
            size_t end = {en};
            values = make_pair(start, end);
            ifsmap[key] = values;
        }
    } else std::cout << "error opening the ifs" << endl;

    ifs.clear();
    ifs.seekg(0, ios::beg);

    double P = TrueOverlapping(ifs); // computing true overlapping reads pairs from fastq
    ifs.close();

    for(size_t i = 0; i < A.cols; ++i) 
    { 
        for(size_t j = A.colptr[i]; j < A.colptr[i+1]; ++j) 
        {
            if(A.values[j] >= S)
            {
                TPandFP++;
                alignment_length = ComputeLength(ifsmap, i, A.rowids[j]); // compute the overlap length between potential overlapping reads pairs
            
                if(alignment_length >= KMER_LENGTH)    
                    TP++;
            }
        }
    }
    cout << "S (minimum number of shared k-mers) = " << S << endl;
    cout << "P = " << P << endl;
    cout << "TP+FP = " << TPandFP << endl;
    cout << "TP = " << TP << endl;
    cout << "Recall = " << TP/P << endl;
    cout << "Precision = " << TP/TPandFP << endl;
    
}
