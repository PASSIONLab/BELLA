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
#include <set>
#include <memory>

using namespace std;

/* Compute the number of true overlapping reads */
double trueOv(ifstream & truth, bool simulated, int ovl_len) 
{
    vector<Interval<std::string>> intervals;
    vector<Interval<std::string>> queries;
    vector<Interval<std::string>>::iterator q;
    double trueoverlaps;

    if(simulated)
    {
        if(truth.is_open()) {
    
            std::string read;
            std::string dontcare2;
            std::string dontcare1;
            std::string ref;
            int start;
            int end;
            while(truth >> ref >> start >> end >> read >> dontcare1 >> dontcare2)
            {
                
                intervals.push_back(Interval<std::string>(start, end, read));
                queries.push_back(Interval<std::string>(start, end, read));
            }
    
        } else std::cout << "Error opening the ground truth file" << endl;
    }
    else
    {
        if(truth.is_open()) {
    
            std::string read;
            int start;
            int end;
            while(truth >> read >> start >> end)
            {
                
                intervals.push_back(Interval<std::string>(start, end, read));
                queries.push_back(Interval<std::string>(start, end, read));
            }
    
        } else std::cout << "Error opening the ground truth file" << endl;
    }

    IntervalTree<std::string> tree;
    vector<size_t> treecounts;

    tree = IntervalTree<std::string>(intervals);

    for (q = queries.begin(); q != queries.end(); ++q) 
    {
        vector<Interval<std::string>> results;
        tree.findOverlapping(q->start, q->stop, results, ovl_len);
        treecounts.push_back(results.size());
    }

    for(size_t t = 0; t < treecounts.size(); ++t) {
        trueoverlaps = trueoverlaps + (double)treecounts[t];
    }

    return trueoverlaps;
}

/* Compute the overlap length between potential overlapping reads pairs */
int computeLength(map<string, pair<int,int>> & sammap, string & col_nametag, string & row_nametag) 
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

void benchmarkingAl(ifstream & groundtruth, ifstream & bella, ifstream & minimap, ifstream & mhap, 
    ifstream & blasr, ifstream & daligner, bool simulated, int ovl_len, int nread) // add blasr && daligner && mhap && bella
{
    map<string, pair<int,int>> ovlsmap;
    map<pair<string,string>, bool> check_bella;
    map<pair<string,string>, bool> check_mini;
    map<pair<string,string>, bool> check_mhap;
    map<pair<string,string>, bool> check_blasr;
    map<pair<string,string>, bool> check_dal;
    map<pair<string,string>, bool>::iterator it;
    
    int alignment_length;
    
    double ovlsbella = 0, truebella = 0;
    double ovlsminimap = 0, trueminimap = 0;
    double ovlsmhap = 0, truemhap = 0;
    double ovlsblasr = 0, trueblasr = 0;
    double ovlsdal = 0, truedal = 0;

    cout << "\nBuilding the ground truth..." << endl;

    if(simulated)
    {
        if(groundtruth.is_open())
        {
            std::string ref;
            std::string read;
            std::string dontcare1;
            std::string dontcare2;
            int start;
            int end;
            std::pair<int,int> coords;
            while(groundtruth >> ref >> start >> end >> read >> dontcare1 >> dontcare2)
            {
                std::string nametag = "@" + read;
                int start_v = start;
                int end_v = end;
                coords = make_pair(start_v, end_v);
                ovlsmap.insert(std::pair<std::string, std::pair<int, int>>(nametag,coords));
            }
        } else std::cout << "Error opening the ground truth file" << endl;
    }
    else
    {
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
        } else std::cout << "Error opening the ground truth file" << endl;
    }
    
    groundtruth.clear();
    groundtruth.seekg(0, ios::beg);

    cout << "Computing BELLA metrics..." << endl;
    if(bella.is_open())
    {
        std::string line;
        while(getline(bella, line))
        {
            ovlsbella++;
            std::stringstream lineStream(line);
            std::string col_nametag, row_nametag;

            getline(lineStream, col_nametag, ' ');
            getline(lineStream, row_nametag, ' ');

            it = check_bella.find(make_pair(col_nametag, row_nametag));
            if(it == check_bella.end())
            {
                check_bella.insert(make_pair(make_pair(col_nametag, row_nametag), true));
                /* Compute the overlap length between potential overlapping reads pairs */
                alignment_length = computeLength(ovlsmap, col_nametag, row_nametag);
                if(alignment_length >= ovl_len)
                    truebella++;
            }

        }
    }

    cout << "Computing Minimap metrics..." << endl;
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

            it = check_mini.find(make_pair(col_nametag, row_nametag));
            if(it == check_mini.end())
            {
                check_mini.insert(make_pair(make_pair(col_nametag, row_nametag), true));
                /* Compute the overlap length between potential overlapping reads pairs */
                alignment_length = computeLength(ovlsmap, col_nametag, row_nametag);
                if(alignment_length >= ovl_len)
                    trueminimap++;
            }
        }
    }

    // as -S option count overlaps only once 
    // (A ov B, but not B ov A), while all other 
    // (included the ground truth) count all ovls and the self-ovls
    double trueminimap_tuned = trueminimap*2+nread;

    cout << "Computing MHAP sensitive metrics..." << endl;
    if(mhap.is_open())
    {
        std::string line;
        while(getline(mhap, line))
        {
            ovlsmhap++;
            std::stringstream lineStream(line);
            std::string col_nametag, row_nametag;

            getline(lineStream, col_nametag, ' ');
            getline(lineStream, row_nametag, ' ');

            col_nametag = "@" + col_nametag;
            row_nametag = "@" + row_nametag;

            it = check_mhap.find(make_pair(col_nametag, row_nametag));
            if(it == check_mhap.end())
            {
                check_mhap.insert(make_pair(make_pair(col_nametag, row_nametag), true));
                /* Compute the overlap length between potential overlapping reads pairs */
                alignment_length = computeLength(ovlsmap, col_nametag, row_nametag);
                if(alignment_length >= ovl_len)
                    truemhap++;
            }
        }
    }

    cout << "Computing BLASR metrics..." << endl;
    if(blasr.is_open())
    {
        std::string line;
        while(getline(blasr, line))
        {
            ovlsblasr++;
            stringstream lineStream(line);
            string col_nametag, row_nametag;

            getline(lineStream, col_nametag, ' ');
            getline(lineStream, row_nametag, ' ');

            col_nametag = "@" + col_nametag;
            row_nametag = "@" + row_nametag;

            it = check_blasr.find(make_pair(col_nametag, row_nametag));
            if(it == check_blasr.end())
            {
                check_blasr.insert(make_pair(make_pair(col_nametag, row_nametag), true));
                /* Compute the overlap length between potential overlapping reads pairs */
                alignment_length = computeLength(ovlsmap, col_nametag, row_nametag);
                if(alignment_length >= ovl_len)
                    trueblasr++;
            }
        }
    }

    cout << "Computing DALIGNER metrics..." << endl;
    if(daligner.is_open())
    {
        std::string line;
        while(getline(daligner, line))
        {
            ovlsdal++;
            stringstream lineStream(line);
            string col_nametag, row_nametag;

            getline(lineStream, col_nametag, ' ');
            getline(lineStream, row_nametag, ' ');

            col_nametag = col_nametag;
            row_nametag = row_nametag;

            it = check_dal.find(make_pair(col_nametag, row_nametag));
            if(it == check_dal.end())
            {
                check_dal.insert(make_pair(make_pair(col_nametag, row_nametag), true));
                /* Compute the overlap length between potential overlapping reads pairs */
                alignment_length = computeLength(ovlsmap, col_nametag, row_nametag);
                if(alignment_length >= ovl_len)
                    truedal++;
            }
        }
    }

    groundtruth.clear();
    groundtruth.seekg(0, ios::beg);
    double truetruth = trueOv(groundtruth, simulated, ovl_len);

    bella.close();
    minimap.close();
    mhap.close();
    blasr.close();
    daligner.close();
    groundtruth.close();

    /* Ground Truth */
    cout << "True overlapping from ground truth = " << truetruth << "\n" << endl;
    /* BELLA Recall and precision */
    cout << "Overlapping from BELLA = " << ovlsbella << endl;
    cout << "True overlapping from BELLA = " << truebella << endl;
    cout << "Recall BELLA = " << truebella/truetruth << endl;
    cout << "Precision BELLA = " << truebella/ovlsbella << "\n" << endl;
    /* Minimap Recall and precision */ 
    cout << "Overlapping from minimap = " << ovlsminimap << endl;
    cout << "True overlapping from minimap = " << trueminimap << endl;
    cout << "Recall minimap = " << trueminimap_tuned/truetruth << endl;
    cout << "Precision minimap = " << trueminimap/ovlsminimap << "\n" << endl;
    /* MHAP Recall and precision */ 
    cout << "Overlapping from MHAP = " << ovlsmhap << endl;
    cout << "True overlapping from MHAP= " << truemhap << endl;
    cout << "Recall MHAP = " << truemhap/truetruth << endl;
    cout << "Precision MHAP = " << truemhap/ovlsmhap << "\n" << endl;
    /* BLASR Recall and precision */ 
    cout << "Overlapping from BLASR = " << ovlsblasr << endl;
    cout << "True overlapping from BLASR = " << trueblasr << endl;
    cout << "Recall BLASR = " << trueblasr/truetruth << endl;
    cout << "Precision BLASR = " << trueblasr/ovlsblasr << "\n" << endl;
    /* DALIGNER Recall and precision */ 
    cout << "Overlapping from DALIGNER = " << ovlsdal << endl;
    cout << "True overlapping from DALIGNER = " << truedal << endl;
    cout << "Recall DALIGNER = " << truedal/truetruth << endl;
    cout << "Precision DALIGNER = " << truedal/ovlsdal << "\n" << endl;
}
