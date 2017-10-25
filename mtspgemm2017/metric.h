#include "CSC.h"
#include "utility.h"
#include "IntervalTree.h"
#include "BitMap.h"
#include "global.h"
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
#include <memory>

#define ERR .15 
#define KMER_LENGTH 17

using namespace std;

/* Compute the number of true overlapping reads */
double trueOv(ifstream & fastafai) 
{
    vector<Interval<std::string>> intervals;
    vector<Interval<std::string>> queries;
    vector<Interval<std::string>>::iterator q;
    double trueoverlaps;

    if(fastafai.is_open()) {

        std::string str;
        int start;
        int end;
        while(fastafai >> str >> start >> end)
        {
            
            intervals.push_back(Interval<std::string>(start, end, str));
            queries.push_back(Interval<std::string>(start, end, str));
        }

    } else std::cout << "Error opening the .sam" << endl;

    IntervalTree<std::string> tree;
    vector<size_t> treecounts;

    tree = IntervalTree<std::string>(intervals);

    for (q = queries.begin(); q != queries.end(); ++q) 
    {
        vector<Interval<std::string>> results;
        tree.findOverlapping(q->start, q->stop, results);
        treecounts.push_back(results.size());
    }

    for(size_t t = 0; t < treecounts.size(); ++t) {
        trueoverlaps = trueoverlaps + (double)treecounts[t];
    }

    return trueoverlaps;
}

/* Compute the overlap length between potential overlapping reads pairs */
int computeLength(map<std::string, std::pair<int,int>> & sammap, std::string & col_nametag, std::string & row_nametag) 
{

    int alignment_length = 0;

    std::map<std::string, std::pair<int,int>>::iterator jit;
    std::map<std::string, std::pair<int,int>>::iterator iit;

    jit = sammap.find(col_nametag); // col name
    iit = sammap.find(row_nametag); // row name 
    
    if(iit->second.first < jit->second.first) {
        if(iit->second.second > jit->second.first) {
            alignment_length = min((iit->second.second - jit->second.first), (jit->second.second - jit->second.first));
        }
    }
    else if (iit->second.first > jit->second.first) {
        if(jit->second.second > iit->second.first) {
            alignment_length = min((jit->second.second - iit->second.first), (iit->second.second - iit->second.first));
        }
    } else alignment_length = min((jit->second.second - iit->second.first), (iit->second.second - iit->second.first)); 

    return alignment_length;
}

void getMetrics(std::ifstream & filename) 
{
    std::ifstream sam("GCF_0001.axt");
    //std::ifstream sam_rv("parsedSamrv.axt");
    std::map<std::string, std::pair<int,int>> sammap;
    //std::map<std::string, std::pair<int,int>> sammaprv;
    int alignment_length;
    double TPandFP = 0, TP = 0;

    if(sam.is_open())
    {
        std::string str;
        int start;
        int end;
        std::pair<int,int> coords;
        while(sam >> str >> start >> end)
        {
            std::string nametag = str;
            int start_v = start;
            int end_v = end;
            coords = make_pair(start_v, end_v);
            sammap.insert(std::pair<std::string, std::pair<int, int>>(nametag,coords));
        }
    } else std::cout << "Error opening the .sam" << endl;

    //if(sam_rv.is_open())
    //{
    //    std::string str;
    //    int start;
    //    int end;
    //    std::pair<int,int> coords;
    //    while(sam_rv >> str >> start >> end)
    //    {
    //        std::string nametag = str;
    //        int start_v = start;
    //        int end_v = end;
    //        coords = make_pair(start_v, end_v);
    //        sammaprv.insert(std::pair<std::string, std::pair<int, int>>(nametag,coords));
    //    }
    //} else std::cout << "Error opening the .sam-rv" << endl;

    sam.clear();
    sam.seekg(0, ios::beg);

    //sam_rv.clear();
    //sam_rv.seekg(0, ios::beg);

    if(filename.is_open())
    {
        std::string line;
        while(getline(filename, line))
        {
            TPandFP++;
            std::stringstream lineStream(line);
            std::string col_nametag, row_nametag;

            getline(lineStream, col_nametag, ' ');
            getline(lineStream, row_nametag, ' ');

            /* Compute the overlap length between potential overlapping reads pairs */
            alignment_length = computeLength(sammap, col_nametag, row_nametag);
            if(alignment_length >= 17)
                TP++;

           // alignment_length = computeLength(sammaprv, col_nametag, row_nametag);
           // if(alignment_length >= 17)
           //     TP++;
        }
    }

    sam.clear();
    sam.seekg(0, ios::beg);

    //sam_rv.clear();
    //sam_rv.seekg(0, ios::beg);

    /* Compute true overlapping reads pairs from fastq */
    //double P_sam = trueOv(sam);
    //double P_sam_rv = trueOv(sam_rv);

    double P = trueOv(sam);

    filename.close();
    sam.close();
    // sam_rv.close();

    cout << "Overlapping from BELLA = " << TPandFP << endl;
    cout << "True overlapping from fastq = " << P << endl;
    cout << "True overlapping from BELLA = " << TP << endl;
    cout << "Recall = " << TP/P << endl;
    cout << "Precision = " << TP/TPandFP << endl; 
}
