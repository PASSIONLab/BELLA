#ifndef _BELLA_MINIMIZER_H_
#define _BELLA_MINIMIZER_H_

#include "../kmercode/hash_funcs.h"
#include "../kmercode/Kmer.hpp"
#include <queue>
#include <vector>
#include <utility>

int robustwinnow = 1;

void furtherPop(std::deque< std::pair<int, uint64_t> > & deq)
{
    if(robustwinnow)
    {
         while (deq.size() > 1 && deq.front().second == deq[1].second)
         {
             deq.pop_front();
         }
    }
}

uint64_t getOrder(const Kmer & mer)
{
    return mer.rep().hash();
}

void sample(std::pair<int, uint64_t> min, std::vector< int > & output)
{
    if(output.empty() || min.first != output.back())
        output.push_back(min.first);
}


void furtherSample(std::deque< std::pair<int, uint64_t> > & deq, std::vector< int > & output)
{
    if(!robustwinnow)
    {
        while (deq.size() > 1 && deq.front().second == deq[1].second)
        {
            deq.pop_front();
            sample(deq.front(), output);
        }
    }
}

// this is a canonical strand minimizer
// to avoid strand ambiguity, we advise choosing an odd number for kmer length
void getMinimizers(size_t window, const std::vector<Kmer>& input, std::vector<int>& output)
{
    std::deque<std::pair<int, uint64_t>> deq; // this double-ended queue will naturally be sorted by construction

    size_t total = input.size();
    for (size_t i = 0; i < total; ++i)
    {

        std::pair<int, uint64_t> qentry = std::make_pair(i, getOrder(input[i]));
        while (!deq.empty() && deq.back().second > qentry.second)
        {
            deq.pop_back();
        }
        deq.push_back(qentry);
	
	while(!deq.empty() && deq.front().first <= static_cast<int>(i)-window)
        {
            furtherPop(deq);  // discard all identical kmers to first
            deq.pop_front();  // discard first out-of-range k-mer
        }

        if (!deq.empty())
        {
            sample(deq.front(), output);
        }

        furtherSample(deq, output);
    }
}

#endif


