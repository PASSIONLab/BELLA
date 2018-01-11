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

#define ERR .15 
#define THR 2000

using namespace std;

/* Compute the number of true overlapping reads */
double trueOv(ifstream & truth) 
{
    vector<Interval<std::string>> intervals;
    vector<Interval<std::string>> queries;
    vector<Interval<std::string>>::iterator q;
    double trueoverlaps;

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

///* Compute the overlap length between potential overlapping reads pairs */
//int countTrue(map<std::string, std::pair<int,int>> & sammap, std::set<pair<std::string, std::string>> & set) 
//{
//    int alignment = 0;
//    int count = 0;
//
//    std::map<std::string, std::pair<int,int>>::iterator jit;
//    std::map<std::string, std::pair<int,int>>::iterator iit;
//
//    for (std::set<pair<std::string, std::string>>::iterator it=set.begin(); it!=set.end(); ++it) {
//
//        jit = sammap.find(it->first);   // col name
//        iit = sammap.find(it->second);  // row name 
//
//        if (jit != sammap.end() && iit != sammap.end()) {
//    
//            if(iit->second.first < jit->second.first) {
//                if(iit->second.second > jit->second.first) {
//                    alignment = min((iit->second.second - jit->second.first), (jit->second.second - jit->second.first));
//                }
//            }
//            else if (iit->second.first > jit->second.first) {
//                if(jit->second.second > iit->second.first) {
//                    alignment = min((jit->second.second - iit->second.first), (iit->second.second - iit->second.first));
//                }
//            } else alignment = min((jit->second.second - iit->second.first), (iit->second.second - iit->second.first)); 
//
//            if(alignment >= THR)
//                count++;
//        }
//    }
//    return count;
//}

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

void benchmarkingAl(std::ifstream & groundtruth, std::ifstream & bella, std::ifstream & minimap, std::ifstream & mhap, std::ifstream & blasr, std::ifstream & daligner) // add blasr && daligner && mhap && bella
{
    std::map<std::string, std::pair<int,int>> ovlsmap;
    std::map<std::pair<std::string,std::string>, bool> check_bella;
    std::map<std::pair<std::string,std::string>, bool> check_mini;
    std::map<std::pair<std::string,std::string>, bool> check_mhap;
    //std::map<std::pair<std::string,std::string>, bool> check_mhap_al;
    std::map<std::pair<std::string,std::string>, bool> check_blasr;
    std::map<std::pair<std::string,std::string>, bool> check_dal;
    std::map<std::pair<std::string,std::string>, bool>::iterator it;
    int alignment_length;
    double ovlsbella = 0, truebella = 0;
    double ovlsminimap = 0, trueminimap = 0;
    double ovlsmhap = 0, truemhap = 0;
    // double ovlsmhap_al = 0, truemhap_al = 0;
    double ovlsblasr = 0, trueblasr = 0;
    double ovlsdal = 0, truedal = 0;

    cout << "\nBuilding the ground truth from BWA-MEM filtered sam..." << endl;

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

    groundtruth.clear();
    groundtruth.seekg(0, ios::beg);

    //int countdown = ovlsmap.size();
    //cout << ovlsmap.size() << endl;
//
    //for(std::map<std::string, std::pair<int,int>>::iterator i=ovlsmap.begin(); i!=ovlsmap.end(); i++) {
    //    for(std::map<std::string, std::pair<int,int>>::iterator j=ovlsmap.begin(); j!=ovlsmap.end(); j++) {
    //        set_bwam.insert(make_pair(i->first, j->first));
    //        set_bwam.insert(make_pair(j->first, i->first));
    //    }
    //    countdown--;
    //    cout << countdown << endl; 
    //}
//
    //int truetruth = countTrue(ovlsmap, set_bwam);

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
                if(alignment_length >= THR)
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
                if(alignment_length >= THR)
                    trueminimap++;
            }
        }
    }

    // as -S option count overlaps only once 
    // (A ov B, but not B ov A), while all other 
    // (included the ground truth) count all ovls and the self-ovls
    double trueminimap_tuned = trueminimap*2 + 16890;

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
                if(alignment_length >= THR)
                    truemhap++;
            }
        }
    }

    /*cout << "Computing MHAP sensitive+alignment metrics..." << endl;

    if(mhap_al.is_open())
    {
        std::string line;
        while(getline(mhap_al, line))
        {
            ovlsmhap_al++;
            std::stringstream lineStream(line);
            std::string col_nametag, row_nametag;

            getline(lineStream, col_nametag, ' ');
            getline(lineStream, row_nametag, ' ');

            col_nametag = "@" + col_nametag;
            row_nametag = "@" + row_nametag;

            it = check_mhap_al.find(make_pair(col_nametag, row_nametag));
            if(it == check_mhap_al.end())
            {
                check_mhap_al.insert(make_pair(make_pair(col_nametag, row_nametag), true));
                // Compute the overlap length between potential overlapping reads pairs 
                alignment_length = computeLength(ovlsmap, col_nametag, row_nametag);
                if(alignment_length >= THR)
                    truemhap_al++;
            }
        }
    }*/

    cout << "Computing BLASR metrics..." << endl;
    if(blasr.is_open())
    {
        std::string line;
        while(getline(blasr, line))
        {
            ovlsblasr++;
            std::stringstream lineStream(line);
            std::string col_nametag, row_nametag;

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
                if(alignment_length >= THR)
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
            std::stringstream lineStream(line);
            std::string col_nametag, row_nametag;

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
                if(alignment_length >= THR)
                    truedal++;
            }
        }
    }

    groundtruth.clear();
    groundtruth.seekg(0, ios::beg);
    double truetruth = trueOv(groundtruth);

    bella.close();
    minimap.close();
    mhap.close();
    blasr.close();
    daligner.close();
    groundtruth.close();

    /* Ground Truth */
    cout << "True overlapping from BWA-MEM = " << truetruth << "\n" << endl;
    /* BELLA Recall and precision */
    cout << "Overlapping from BELLA = " << ovlsbella << endl;
    cout << "True overlapping from BELLA = " << truebella << endl;
    cout << "Recall BELLA = " << truebella/(double)truetruth << endl;
    cout << "Precision BELLA = " << truebella/ovlsbella << "\n" << endl;
    /* Minimap Recall and precision */ 
    cout << "Overlapping from minimap = " << ovlsminimap << endl;
    cout << "True overlapping from minimap = " << trueminimap << endl;
    cout << "Recall minimap = " << trueminimap_tuned/(double)truetruth << endl;
    cout << "Precision minimap = " << trueminimap/ovlsminimap << "\n" << endl;
    /* MHAP Recall and precision */ 
    cout << "Overlapping from MHAP = " << ovlsmhap << endl;
    cout << "True overlapping from MHAP= " << truemhap << endl;
    cout << "Recall MHAP = " << truemhap/(double)truetruth << endl;
    cout << "Precision MHAP = " << truemhap/ovlsmhap << "\n" << endl;
    /* MHAP+alignement Recall and precision */ 
    //cout << "Overlapping from MHAP+alignement = " << ovlsmhap_al << endl;
    //cout << "True overlapping from MHAP+alignement = " << truemhap_al << endl;
    //cout << "Recall MHAP+alignement = " << truemhap_al/(double)truetruth << endl;
    //cout << "Precision MHAP+alignement = " << truemhap_al/ovlsmhap_al << "\n" << endl;
    /* BLASR Recall and precision */ 
    cout << "Overlapping from BLASR = " << ovlsblasr << endl;
    cout << "True overlapping from BLASR = " << trueblasr << endl;
    cout << "Recall BLASR = " << trueblasr/truetruth << endl;
    cout << "Precision BLASR = " << trueblasr/ovlsblasr << "\n" << endl;
    /* DALIGNER Recall and precision */ 
    cout << "Overlapping from DALIGNER = " << ovlsdal << endl;
    cout << "True overlapping from DALIGNER = " << truedal << endl;
    cout << "Recall DALIGNER = " << truedal/truetruth << endl;
    cout << "Precision DALIGNER = " << truedal/ovlsdal << endl;
}
