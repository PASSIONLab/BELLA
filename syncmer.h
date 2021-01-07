#ifndef _BELLA_SYNCMER_H_
#define _BELLA_SYNCMER_H_

#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include <queue>
#include <vector>
#include <utility>



// int robustwinnow = 1;

// void furtherPop(std::deque< std::pair<int, uint64_t> > & deq)
// {
//     if(robustwinnow)
//     {
//          while (deq.size() > 1 && deq.front().second == deq[1].second)
//          {
//              deq.pop_front();
//          }
//     }
// }

// uint64_t getOrder(const Kmer & mer)
// {
//     return mer.rep().hash();
// }

// void sample(std::pair<int, uint64_t> min, std::vector< int > & output)
// {
//     if(output.empty() || min.first != output.back())
//         output.push_back(min.first);
// }


// void furtherSample(std::deque< std::pair<int, uint64_t> > & deq, std::vector< int > & output)
// {
//     if(!robustwinnow)
//     {
//         while (deq.size() > 1 && deq.front().second == deq[1].second)
//         {
//             deq.pop_front();
//             sample(deq.front(), output);
//         }
//     }
// }

const int smerlen = 5; //make it input param

bool isSyncmer(const Kmer & mer, int &kmerlen)
{
    string kmer = mer.toString();
    Kmer stKmer(kmer.substr(0, smerlen).c_str(), smerlen); // first s-mer of the k-mer
    Kmer endKmer(kmer.substr(kmerlen - smerlen, smerlen).c_str(), smerlen); // last s-mer of the k-mer

    uint64_t stsmer =  stKmer.hash(); //(kmer >> (2*(31-(smerlen - 1)))) & mask; // first s-mer of the k-mer
    uint64_t endsmer = endKmer.hash();    
    
    for (int i = 1; i < kmerlen - smerlen ; ++i)
    {
        Kmer tmpsmer(kmer.substr(i, kmerlen).c_str(), smerlen);

        if(tmpsmer.hash() < stsmer && tmpsmer.hash() < endsmer){
            return false;
        }
    }           
    return true;

}

// this is a canonical strand minimizer
// to avoid strand ambiguity, we advise choosing an odd number for kmer length
void getSyncmers(int klen, const std::vector<Kmer> & input, std::vector< int> & output, int ishopc)
{ 
    for (size_t i=0; i< input.size(); ++i)
    {
        if(isSyncmer(input[i], klen))
            output.push_back(i);
    }
    // cout << "#kmer" << total << "#syncmers: " << output.size() << endl;
    // for(auto minmer:output)
    // {
    //    cout << input[minmer] << " ";
    // }
    // cout << endl;
}


#endif


