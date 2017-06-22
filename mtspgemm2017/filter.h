#include "CSC.h"
#include "utility.h"
#include "BitMap.h"
#include "global.h"

extern "C" {
#include "../DALIGNER/align.h"
}

#include <cstdio> 
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
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h> 
#include "../edlib/edlib/include/edlib.h"
#include <sys/types.h> 
#include <sys/stat.h> 
#include <math.h>
#include <limits.h>

#define UPPER_BOUND 7
#define ERR .15 
#define P_CORR_BOTH 0.7225
#define CORR_FACT 1.5

/* Estimation of overlapping region length */
template <typename IT>
int FindOverlap(std::string & aseq, std::string & bseq, IT & apos, IT & bpos, int & alen, int & blen, int & left, int & right) 
{ 
    int atemp, btemp;
    
    if(apos <= bpos)
      left = (int)apos;
    else left = (int)bpos;

    atemp = alen - (int)apos;
    btemp = blen - (int)bpos; 

    if(atemp <= btemp)
        right = atemp;
    else right = btemp;

    int overlap = (left+right)*CORR_FACT;

    if(apos <= bpos){
      aseq = aseq.substr(0, overlap);
      bseq = bseq.substr(bpos-left, overlap);
    }
    else {
      bseq = bseq.substr(0, overlap);
      aseq = aseq.substr(apos-left, overlap);
    }

    return overlap;
}

bool Filter(std::string & aseq, std::string & bseq, int & alen, int & blen) 
{
  bool align;
  EdlibAlignResult result;
  char *c_aseq = new char[alen+1];
  strcpy(c_aseq, aseq.c_str());
  char *c_bseq = new char[blen+1];
  strcpy(c_bseq, bseq.c_str());
  int shr;

  if(alen <= blen)
    shr = alen;
  else shr = blen;

  int ed = shr - (shr*P_CORR_BOTH - shr*P_CORR_BOTH*(1 - P_CORR_BOTH)); // binomial distribution: mean + variance = upper bound for edit distance 
  
  result = edlibAlign(c_aseq, alen, c_bseq, blen, edlibNewAlignConfig(ed, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE));
  delete [] c_aseq;
  delete [] c_bseq;

  if(result.editDistance == -1)
    align = false;
  else align = true;

  edlibFreeAlignResult(result);

  return align;
}

template <typename IT, typename NT>
void LocalAlignment(const CSC<IT,NT> & A, std::vector<string> & reads) 
{
  std::string aseq, bseq;
  int alen, blen, left, right;
  bool align;


	for(size_t i = 0; i < A.cols; ++i) 
  { 
    for(size_t j = A.colptr[i]; j < A.colptr[i+1]; ++j) 
    { 
      if(A.values[j].first < UPPER_BOUND) 	// Pairs sharing more than upper_bound k-mers have high precision
      {
        aseq = reads[A.rowids[j]];
        bseq = reads[i];
        alen = (int)aseq.length();
        blen = (int)bseq.length();
        
        //overlap = FindOverlap(aseq, bseq, A.values[j].second.first, A.values[j].second.second, alen, blen, left, right); // apos, bpos, alen, blen
        align = Filter(aseq, bseq, alen, blen);
        
        if(align == false)
        {
          A.values[j].first = 0;
        }
      }     
    }
  }
} 