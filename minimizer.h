#ifndef _BELLA_MINIMIZER_H_
#define _BELLA_MINIMIZER_H_

#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include <queue>
#include <vector>
#include <utility>



int robustwinnow = 0;

void furtherPop(std::deque< std::pair<int, uint64_t> > & deq)
{
    if(robustwinnow)
    {
         while (deq.front().second == deq[1].second)
         {
             deq.pop_front();
         }
    }
}

uint64_t getOrder(const Kmer & mer)
{
    return mer.rep().hash();
}

void sample(std::pair<int, uint64_t> min, std::vector< std::pair<int, uint64_t> > & output)
{
    if(output.empty() || output.back() != min)
    {
        output.push_back(min);
    }
}


void furtherSample(std::deque< std::pair<int, uint64_t> > & deq, std::vector< std::pair<int, uint64_t> > & output)
{
    if(!robustwinnow)
    {
        while (deq.front().second == deq[1].second)
        {
            sample(deq.front(), output);
            deq.pop_front();
        }
    }
}

// this is a canonical strand minimizer
// to avoid strand ambiguity, we advise choosing an odd number for kmer length
void getMinimizers(int window, const std::vector<Kmer> & input, std::vector< std::pair<int, uint64_t> > & output, int ishopc)
{
    std::deque< std::pair<int, uint64_t> > deq; // this double-ended queue will naturally be sorted by construction

    size_t total = input.size();
    for (size_t i=0; i< total; ++i)
    {
        std::pair<int, uint64_t> qentry = std::make_pair(i, getOrder(input[i]));
        while (!deq.empty() && deq.back().second > qentry.second)
        {
            deq.pop_back();
        }
        deq.push_back(qentry);
        if(deq.front().first <= static_cast<int>(i)-window)
        {
            while(deq.front().first <= static_cast<int>(i)-window)
            {
                deq.pop_front();  // discard out-of-range k-mer
            }
            furtherPop(deq);
        }
        sample(deq.front(), output);

        furtherSample(deq, output);
    }
}


#endif


