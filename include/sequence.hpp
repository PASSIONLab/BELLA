#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include "common/common.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "../kmercode/Kmer.hpp"

using namespace std;

char
complementbase(char n) {
	switch(n)
	{
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	}
	assert(false);
	return ' ';
}

std::string
reversecomplement(const std::string& seq) {

	std::string cpyseq = seq;
	std::reverse(cpyseq.begin(), cpyseq.end());

	std::transform(
		std::begin(cpyseq),
		std::end  (cpyseq),
		std::begin(cpyseq),
	complementbase);

	return cpyseq;
}

// GGGG: function imported from https://github.com/marbl/Winnowmap/blob/391b045d139a97db577b53216905f237a9c33b4c/ext/meryl/src/utility/src/utility/sequence.C

//  Compress homopolymer runs to a single letter.  Returns the length of the
//  compressed sequence.
//
//  'bases' does not need to be NUL terminated.
//
//  If 'compr' is supplied, the compressed sequence is returned.  'compr' can be
//  the same array as 'bases'.  The output string is always NUL terminated.
//
//  If 'ntoc' is supplied, a mapping from 'normal' position to 'compressed'
//  position is returned.  This output is INCORRECT if 'skip' is set (the
//  values of any skipped bases are undefined).
//
//  If 'skip' is set, ignore these bases at the start of the sequence that
//  are the same.  This is to allow one to homopoly compress a sequence in
//  chunks; compress the first 1000 bases, the the next 1000, and so on.
//  Each pass sets skip to the last base of the previous chunk.
//
uint32
homopolyCompress(char *bases, uint32 basesLen, char *compr, uint32 *ntoc, char skip)
{
  uint32  cc = 0;  //  position of the start of the run
  uint32  rr = 1;  //  position of the scan head
  uint32  sl = 0;  //  length of the compressed sequence

  while ((bases[cc] == skip) &&   //  If 'skip' is set, ignore these bases
         (cc < basesLen)) {        //  at the start of 'bases'.
    cc++;
    rr++;
  }

  if (compr)                      //  Either the first base, or
    compr[sl] = bases[cc];        //  the terminating NUL.

  if (ntoc)                       //  Save the mapping from the first
    ntoc[cc] = sl;                //  normal to first compressed base.

  if (basesLen == 0)
    return(0);

  sl++;

  while (rr < basesLen) {

    //  In a run, move the scan head one position, and set the
    //  mapping to the last compressed base.
    if ((bases[cc] | 0x20) == (bases[rr] | 0x20)) {    //  Convert to lowercase before comparing.
      if (ntoc)
        ntoc[rr] = sl - 1;
      rr++;
    }

    //  Just ended a run. Set the start of the (next) run
    //  to the current position, move the current position
    //  to the next base, and increase the length of the
    //  compressed sequence.
    else {
      if (compr)
        compr[sl] = bases[rr];
      if (ntoc)
        ntoc[rr] = sl;
      cc = rr;
      rr++;
      sl++;
    }
  }

  //  Terminate the compressed string.
  if (compr)
    compr[sl] = 0;

  //  The 'space' after the end of the bases maps to the 'space'
  //  after the compressed bases.
  if (ntoc)
    ntoc[rr]  = sl;

  return(sl);
}

#endif // _SEQUENCE_H_