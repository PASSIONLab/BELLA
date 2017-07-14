#include "CSC.h"
#include "utility.h"
#include "IntervalTree.h"
#include "BitMap.h"
#include "global.h"
#include "../edlib/edlib/include/edlib.h"
#include <omp.h>
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
#include <unistd.h>
#include <ctype.h> 

#define ERR .15 
#define KMER_LENGTH 17
using namespace std;

#ifdef _MULTPTR
/* Compute the distance between two k-mers on the same read */
template <typename IT>
double ComputeDistance(IT & fst, IT & snd) 
{

    double delta;
    double min, max, diff;

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

/* Find when a k-mers pair represents an evidence of potential overlap or not */
bool PotentialOverlap(double & dfst, double & dsnd) 
{

    double dmin, dmax; /* Define the shortest and longest distance */
    double Lmin, Lmax; /* Compute the shortest length of the longest distance and the longest length of the shortest distance */
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
    
    Lmax = dmin/p_del; /* Maximum length of the shortest delta given the probability of deletion */
    Lmin = dmax/p_ins; /* Minimum length of the longest delta given the probability of insertion */
    
    if(Lmax > Lmin) {
        region = true;
    } else region = false;

    return region;
}

/* Estimation of overlapping region length
   Reads length has to be included and position of the same k-mer on different reads */
template <typename IT>
double ExpectedKmers(IT & i, IT & j, IT & len_i, IT & len_j) 
{

    double left;    /* Define the left margin */
    double right;   /* Define the right margin */
    double estime, p;  /* Expected overlap region length */
    size_t temp_i, temp_j;
    
    if(i <= j)
        left = (double)i;
    else left = (double)j;

    temp_i = len_i - i;
    temp_j = len_j - j; 

    if(temp_i <= temp_j)
        right = (double)temp_i;
    else right = (double)temp_j;

    estime = left + right; /* Estimated overlap */
    p = 1-pow((1-pow((1-ERR), (2*KMER_LENGTH))), estime); /* Expected number of k-mers */

    return p;
}

/* Reads length has to be included and position of the same k-mer on different reads */
template <typename IT>
double MaxGap(IT & i, IT & j, IT & len_i, IT & len_j) 
{

    double left;    /* Define the left margin */
    double right;   /* Define the right margin */
    double estime;  /* Expected overlap region length */
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

    estime = left + right; /* Estimated overlap */

    return estime;
}
#endif

#ifdef _LOCALIGN
/* EDLIB Local Alignment */
bool Filter(std::string & aseq, std::string & bseq, int & alen, int & blen) 
{
  bool align;
  EdlibAlignResult result;
  char *c_aseq = new char[alen+1];
  strcpy(c_aseq, aseq.c_str());
  char *c_bseq = new char[blen+1];
  strcpy(c_bseq, bseq.c_str());

  result = edlibAlign(c_aseq, alen, c_bseq, blen, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH));
  delete [] c_aseq;
  delete [] c_bseq;

  if(result.alignmentLength < 300)
    align = false;
  else align = true;

  edlibFreeAlignResult(result);
  return align;
}
#endif

/* Compute the number of true overlapping reads */
double TrueOverlapping(ifstream & ifs) 
{

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

/* Compute the overlap length between potential overlapping reads pairs */
int ComputeLength(map<int, pair<int,int>> & ifsmap, int & col, int & row) 
{

    int alignment_length = 0;
    std::map<int, pair<int,int>>::iterator jit;
    std::map<int, pair<int,int>>::iterator iit;

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

void LocalAlignmentTest(std::ifstream & filename) 
{
    std::ifstream ifs("test_01.axt"); /* To be generalized */
    std::map<int, pair<int,int>> ifsmap;
    double TPandFP = 0, TP = 0;
    int alignment_length;

    #ifdef _MULTPTR
    double di; /* k-mers distances on read i */
    double dj; /* k-mers distances on read j */
    #endif

    #ifdef _LOCALIGN
    std::string aseq, bseq;
    int alen, blen;
    double LA = 0;
    #endif

    /* Create reads map from axt file to be used to find potential overlapping reads pairs */
    if(ifs.is_open()) {
        int id;
        int st;
        int en;
        std::pair<int, int> values;
        while(ifs >> id >> st >> en) {
            int key = {id};
            int start = {st};
            int end = {en};
            values = make_pair(start, end);
            ifsmap[key] = values;
        }
    } else std::cout << "error opening the ifs" << endl;

    ifs.clear();
    ifs.seekg(0, ios::beg);

    /* Compute true overlapping reads pairs from fastq */
    double P = TrueOverlapping(ifs);
    ifs.clear();
    ifs.seekg(0, ios::beg);

    if(filename.is_open())
    {
        std::string line;
        while(getline(filename, line))
        {
            TPandFP++;
            std::stringstream lineStream(line);
            std::string col, row;

            getline(lineStream, col, ',');
            getline(lineStream, row, ',');

            int colid = stoi(col);
            int rowid = stoi(row);

            #ifdef _LOCALIGN
            aseq = reads[A.rowids[j]];
            bseq = reads[i];
            alen = (int)aseq.length();
            blen = (int)bseq.length();
            #endif

            /* Compute the overlap length between potential overlapping reads pairs */
            alignment_length = ComputeLength(ifsmap, colid, rowid);
            if(alignment_length >= KMER_LENGTH)
            {
                TP++;

                #ifdef _LOCALIGN
                align = Filter(aseq, bseq, alen, blen);
                if(align == true)
                    LA++;
                #endif
            }
        }
    }
    filename.close();
    ifs.close();

    cout << "TP+FP after multiply = " << TPandFP << endl;
    cout << "P from fastq = " << P << endl;
    cout << "TP from fastq = " << TP << endl;
    cout << "Recall = " << TP/P << endl;
    cout << "Precision = " << TP/TPandFP << endl; 
}
