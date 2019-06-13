#include "IntervalTree.h"
#include "common.h"
//#include <omp.h>
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
#include <sstream>
#include <set>
#include <memory>
#include <typeinfo>

using namespace std;

uint32_t trueOv(truth_m & truthInfo, bool sim, int minOv) 
{
    intervals_v intervals;
    intervals_v queries;
    intervals_v::iterator q;
    truth_m::iterator key; // outer iterator
    info_v::iterator it; // inner iterator
    uint32_t thOv = 0;

    for(key = truthInfo.begin(); key != truthInfo.end(); ++key)
    {
         for(it = key->second.begin(); it != key->second.end(); ++it)
         {
             intervals.push_back(Interval<string>(it->start, it->end, it->read));
             queries.push_back(Interval<string>(it->start, it->end, it->read));
         }

         IntervalTree<string> tree;
         vector<size_t> treeCount;
         
         tree = IntervalTree<string>(intervals); // reference 
         
         for (q = queries.begin(); q != queries.end(); ++q) // tree search for a given reference
         {
             vector<Interval<string>> results;
             tree.findOverlapping(q->start, q->stop, q->value, results, minOv);
             treeCount.push_back(results.size());
         }
         
         for(size_t t = 0; t < treeCount.size(); ++t) // count for a given reference
         { 
             thOv = thOv + (double)treeCount[t];  // cumulative
         }

         intervals.clear();
         queries.clear();
    }
    return thOv;
}

int64_t getTruth(ifstream & th, bool sim, int minOv, mmap_m & seqmap, uint32_t numth)
{
    truth_m isInThere;
    truth_m::iterator iter;
    mecat_m namestable;

    cout << "\nBuilding the ground truth ..." << endl;

    if(sim) // simulated data
        if(th.is_open())
        {
            string ref;
            string prev;
            string read;
            string skip1;
            string skip2;
            refInfo ovlInfo;
            readInfo perRead;
            int start;
            int end;

            while(th >> ref >> start >> end >> read >> skip1 >> skip2)
            {
                perRead.ref = ref;
                perRead.start = start;
                perRead.end = end;
                seqmap.insert(make_pair("@"+read, perRead));

                ovlInfo.ref = ref;
                ovlInfo.read = "@" + read;
                ovlInfo.start = start;
                ovlInfo.end = end;

                iter = isInThere.find(ref);
                if(iter == isInThere.end())
                {
                    info_v temp;
                    temp.push_back(ovlInfo); // all the element of a chromosome
                    isInThere.insert(truth_m::value_type(ref,temp));
                }
                else
                {
                    iter->second.push_back(ovlInfo);
                    isInThere[ref] = iter->second;
                }
            }
            cout << "numReads: " << seqmap.size() << endl;
            cout << "numSubreferences: " << isInThere.size() << endl;
        } else cout << "Error opening the ground truth" << endl;
    else // real data
        if(th.is_open())
        {
            refInfo ovlInfo;
            readInfo perRead;
            string prev;
            string read;
            string ref;
            int start;
            int end;

            while(th >> ref >> read >> start >> end)
            {
                perRead.ref = ref;
                perRead.start = start;
                perRead.end = end;
                seqmap.insert(make_pair("@"+read,perRead));

                ovlInfo.ref = ref;
                ovlInfo.read = "@" + read;
                ovlInfo.start = start;
                ovlInfo.end = end;

                iter = isInThere.find(ref);
                if(iter == isInThere.end())
                {
                    info_v temp;
                    temp.push_back(ovlInfo); // all the element of a chromosome
                    isInThere.insert(truth_m::value_type(ref,temp));
                }
                else
                {
                    iter->second.push_back(ovlInfo);
                    isInThere[ref] = iter->second;
                }
            }
            cout << "numReads: " << seqmap.size() << endl;
            cout << "numSubreferences: " << isInThere.size() << endl;
        } else cout << "Error opening the ground truth" << endl;

    th.clear();
    th.seekg(0,ios::beg);

    numth = trueOv(isInThere, sim, minOv);
    cout << "\nPairs in ground truth > " << minOv << "bp = " << numth << "\n" << endl;
    return numth;
}
