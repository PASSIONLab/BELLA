#ifndef _BELLA_MINIMIZER_H_
#define _BELLA_MINIMIZER_H_

#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include <queue>


// this is a canonical strand minimizer
// to avoid strand ambiguity, we advise choosing an odd number for kmer length
void GetMinimizers(int window, const vector<Kmer> & input, vector<Kmer> & output, int ishopc)
{
    std::deque<int> deq;

    size_t total = input.size();
    for (size_t i=0; i< total; ++i)
    {
        pair<size_t, uint64_t> key = make_pair(i, input[i].rep().hash());
        while (!deq.empty() && deq.back().second > key.second)
        {
            deq.pop_back();
        }
    }
}



#endif


