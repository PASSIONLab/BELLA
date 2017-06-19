#include "CSC.h"
#include "utility.h"
#include "align.h"
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
#include <bitset>
#include <map>

#define UPPER_BOUND 5

// Estimation of overlapping region length
template <typename IT>
void FindLowAndHigh(IT & apos, IT & bpos, IT & alen, IT & blen, IT & left, IT & right) 
{ 
    double estime;  // expected overlap region length
    size_t atemp, btemp;
    
    if(apos <= bpos)
        left = (double)apos;
    else left = (double)bpos;

    atemp = alen - apos;
    btemp = blen - bpos; 

    if(atemp <= btemp)
        right = (double)atemp;
    else right = (double)btemp;
}

/* TODO: Implement local alignment via bitvector */
bool LAfilter() {}

template <typename IT, typename NT>
void LocalAlignment(const CSC<IT,NT> & A, std::vector<string> & seqs) 
{
  std::string aseq, bseq;
  size_t alen, blen, left, right;
  bool keep;

	for(size_t i = 0; i < A.cols; ++i) 
  { 
    for(size_t j = A.colptr[i]; j < A.colptr[i+1]; ++j) 
    {
      if(A.values[j].first < UPPER_BOUND) 	// Pairs sharing more than upper_bound k-mers have high precision
      {
        aseq = seqs[A.rowids[j]];
        bseq = seqs[i];
        alen = aseq.length();
        blen = blen.length();

        FindLowAndHigh(A.values[j].second.first, A.values[j].second.second, alen, blen, left, right); // apos, bpos, alen, blen
        keep = LAfilter(aseq, bseq, alen, blen, left, right);

        if(keep == false)
          // remove reads pairs, free A.values[j]
      }
    }
  }
}