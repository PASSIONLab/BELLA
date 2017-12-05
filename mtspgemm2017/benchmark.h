#include "IntervalTree.h"
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
#include <sstream>
#include <memory>

#define ERR .15 
#define THR 2000
#define KMER_LENGTH 17

using namespace std;

/* Compute the number of true overlapping reads */
double trueOv(ifstream & truth) 
{
    vector<Interval<std::string>> intervals;
    vector<Interval<std::string>> queries;
    vector<Interval<std::string>>::iterator q;
    double trueoverlaps;

    if(truth.is_open()) {

        //std::string ref;
        std::string read;
        //std::string dontcare1;
        //std::string dontcare2;
        //std::string dontcare3;
        int start;
        int end;
        while(truth >> read >> start >> end)
        {
            
            intervals.push_back(Interval<std::string>(start, end, read));
            queries.push_back(Interval<std::string>(start, end, read));
        }

    } else std::cout << "Error opening the .axt" << endl;

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

void benchmarkingAl(std::ifstream & groundtruth, std::ifstream & minimap) // add blasr && daligner && mhap && bella
{
    std::map<std::string, std::pair<int,int>> ovlsmap;
    int alignment_length;
    // double ovlsbella = 0, truebella = 0;
    double ovlsminimap = 0, trueminimap = 0;
    // double ovlsblasr = 0, trueblasr = 0;
    // double ovlsdaligner = 0, truedaligner = 0;
    // double ovlsmhap = 0, truemahpa = 0;

    if(groundtruth.is_open())
    {
        std::string read;
        int start;
        int end;
        std::pair<int,int> coords;
        while(groundtruth >> read >> start >> end)
        {
            std::string nametag = "@" + read;
            int start_v = start;
            int end_v = end;
            coords = make_pair(start_v, end_v);
            ovlsmap.insert(std::pair<std::string, std::pair<int, int>>(nametag,coords));
        }
    } else std::cout << "Error opening the .axt" << endl;

    groundtruth.clear();
    groundtruth.seekg(0, ios::beg);

    // if(bella.is_open())
    // {
    //     std::string line;
    //     while(getline(bella, line))
    //     {
    //         ovlsbella++;
    //         std::stringstream lineStream(line);
    //         std::string col_nametag, row_nametag;

    //         getline(lineStream, col_nametag, ' ');
    //         getline(lineStream, row_nametag, ' ');

    //         /* Compute the overlap length between potential overlapping reads pairs */
    //         alignment_length = computeLength(ovlsmap, col_nametag, row_nametag);
    //         if(alignment_length >= THR)
    //             truebella++;
    //     }
    // }

    if(minimap.is_open())
    {
        std::string line;
        while(getline(minimap, line))
        {
            ovlsminimap++;
            std::stringstream lineStream(line);
            std::string col_nametag, row_nametag, dontcare1, dontcare2, dontcare3, dontcare4;
    
            getline(lineStream, col_nametag, '\t' );
            getline(lineStream, dontcare1, '\t' );
            getline(lineStream, dontcare2, '\t' );
            getline(lineStream, dontcare3, '\t' );
            getline(lineStream, dontcare4, '\t' );
            getline(lineStream, row_nametag, '\t' );
    
            col_nametag = "@" + col_nametag;
            row_nametag = "@" + row_nametag;
    
            /* Compute the overlap length between potential overlapping reads pairs */
            alignment_length = computeLength(ovlsmap, col_nametag, row_nametag);
            if(alignment_length >= THR)
                trueminimap++;
        }
    }

    //if(mhap.is_open())
    //{
    //    std::string line;
    //    while(getline(mhap, line))
    //    {
    //        ovlsmhap++;
    //        std::stringstream lineStream(line);
    //        std::string col_nametag, row_nametag;

    //        getline(lineStream, col_nametag, ' ');
    //        getline(lineStream, row_nametag, ' ');

    //        col_nametag = "@" + col_nametag;
    //        row_nametag = "@" + row_nametag;

    //        /* Compute the overlap length between potential overlapping reads pairs */
    //        alignment_length = computeLength(ovlsmap, col_nametag, row_nametag);
    //        if(alignment_length >= 1000)
    //            truemhap++;
    //    }
    //}
    
    //if(blasr.is_open())
    //{
    //    std::string line;
    //    while(getline(blasr, line))
    //    {
    //        ovlsblasr++;
    //        std::stringstream lineStream(line);
    //        std::string col_nametag, row_nametag;

    //        getline(lineStream, col_nametag, ' ');
    //        getline(lineStream, row_nametag, ' ');

    //        col_nametag = "@" + col_nametag;
    //        row_nametag = "@" + row_nametag;

    //        /* Compute the overlap length between potential overlapping reads pairs */
    //        alignment_length = computeLength(ovlsmap, col_nametag, row_nametag);
    //        if(alignment_length >= 1000)
    //            trueblasr++;
    //    }
    //} 

    groundtruth.clear();
    groundtruth.seekg(0, ios::beg);
    double truetruth = trueOv(groundtruth);

    // bella.close();
    minimap.close();
    // mhap.close();
    // blasr.close();
    // daligner.close();
    groundtruth.close();

    /* Ground Truth */
    cout << "True overlapping from fastq = " << truetruth << endl;
    /* BELLA Recall and precision */
    //cout << "Overlapping from BELLA = " << ovlsbella << endl;
    //cout << "True overlapping from BELLA = " << truebella << endl;
    //cout << "Recall BELLA = " << truebella/truetruth << endl;
    //cout << "Precision BELLA = " << truebella/ovlsbella << endl;
    /* Minimap Recall and precision */ 
    cout << "Overlapping from minimap = " << ovlsminimap << endl;
    cout << "True overlapping from minimap = " << trueminimap << endl;
    cout << "Recall minimap = " << trueminimap/truetruth << endl;
    cout << "Precision minimap = " << trueminimap/ovlsminimap << endl;
    /* MHAP Recall and precision */ 
    //cout << "Overlapping from BLASR = " << ovlsmhap << endl;
    //cout << "True overlapping from BLASR = " << truemhap << endl;
    //cout << "Recall BLASR = " << truemhap/truetruth << endl;
    //cout << "Precision BLASR = " << truemhap/ovlsmhap << endl;
    /* BLASR Recall and precision */ 
    //cout << "Overlapping from BLASR = " << ovlsblasr << endl;
    //cout << "True overlapping from BLASR = " << trueblasr << endl;
    //cout << "Recall BLASR = " << trueblasr/truetruth << endl;
    //cout << "Precision BLASR = " << trueblasr/ovlsblasr << endl;
    /* DALIGNER Recall and precision */ 
    //cout << "Overlapping from BLASR = " << ovlsdaligner << endl;
    //cout << "True overlapping from BLASR = " << truedaligner << endl;
    //cout << "Recall BLASR = " << truedaligner/truetruth << endl;
    //cout << "Precision BLASR = " << truedaligner/ovlsdaligner << endl;
}
