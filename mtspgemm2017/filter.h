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
#include <sys/types.h> 
#include <sys/stat.h> 
#include <math.h>
#include <limits.h>

#define UPPER_BOUND 5
#define AVE_ERROR .70
#define SPACING 100
#define ERR .15 

/* Estimation of overlapping region length */
template <typename IT>
void FindLowAndHigh(IT & apos, IT & bpos, int & alen, int & blen, int & left, int & right) 
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
}

/* TODO: Implement local alignment via bitvector */
bool Filter(std::string & aseq, std::string & bseq, size_t & apos, size_t & bpos, int & alen, int & blen, int & low, int & hgh) 
{
  Path *path = new Path();
  Alignment *align = new Alignment();
  Work_Data *ework = New_Work_Data();
  Align_Spec *espec = New_Align_Spec(AVE_ERROR, SPACING);
  Path *filteredpath; 

  int anti = (int)apos + (int)bpos;
  int diag = (int)apos - (int)bpos;
  int lbord = diag - diag*ERR;  // issue here
  int hbord = diag + diag*ERR;  // issue here -- seg fault in forward wave after three iteration

  char *caseq = new char[aseq.length()+1];
  strcpy(caseq, aseq.c_str());
  align->aseq = caseq;
  delete [] caseq;

  char *cbseq = new char[bseq.length()+1];
  strcpy(cbseq, bseq.c_str());
  align->bseq = cbseq;
  delete [] cbseq;

  align->alen = alen;
  align->blen = blen;
  align->path = path;
  
  cout << "Before Local_Alignment" << endl;
  filteredpath = Local_Alignment(align, ework, espec, low, hgh, anti, lbord, hbord);
  if(path != NULL)
  {
    cout << "path created" << endl;
    return true;
  } else return false;
}

template <typename IT, typename NT>
void LocalAlignment(const CSC<IT,NT> & A, std::vector<string> & reads) 
{
  std::string aseq, bseq;
  int alen, blen, left = 0, right = 0;
  bool keep = false;

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
        cout << "pos on read i: " << A.values[j].second.first << endl;
        cout << "pos on read j: " << A.values[j].second.second << endl;
        FindLowAndHigh(A.values[j].second.first, A.values[j].second.second, alen, blen, left, right); // apos, bpos, alen, blen
        keep = Filter(aseq, bseq, A.values[j].second.first, A.values[j].second.second, alen, blen, left, right);

        if(keep == false)
          cout << "Remove reads pair" << endl;
      }
    }
  }
}