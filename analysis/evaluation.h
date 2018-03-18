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
#include <unordered_map>
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

struct refInfo {

    string ref;
    string read;
    int start;
    int end;
};

struct readInfo {

    string ref;
    int start;
    int end;
};

typedef vector<refInfo> vectInfo;

/* Compute the number of true overlapping reads */
double trueOv(vector<vectInfo> & truthInfo, bool simulated, int minOvl) 
{
    vector<Interval<std::string>> intervals;
    vector<Interval<std::string>> queries;
    vector<Interval<std::string>>::iterator q;
    vector<vectInfo>::iterator oit; // outer iterator
    vectInfo::iterator iit; // inner iterator
    double trueOvls = 0;

    if(simulated)
    {
       for(oit=truthInfo.begin(); oit!=truthInfo.end(); ++oit)
       {
            for(iit=oit->begin(); iit!=oit->end(); ++iit)
            {
                intervals.push_back(Interval<string>(iit->start, iit->end, iit->read));
                queries.push_back(Interval<string>(iit->start, iit->end, iit->read));
            }

            IntervalTree<string> tree;
            vector<size_t> treeCount;
            
            tree = IntervalTree<string>(intervals); // reference 
            
            for (q = queries.begin(); q != queries.end(); ++q) // tree search for a given reference
            {
                vector<Interval<string>> results;
                tree.findOverlapping(q->start, q->stop, q->value, results, minOvl);
                treeCount.push_back(results.size());
            }
            
            for(size_t t = 0; t < treeCount.size(); ++t) // count for a given reference
            { 
                trueOvls = trueOvls + (double)treeCount[t];  // cumulative
            }

            intervals.clear();
            queries.clear();
       }          
    }
    //else // need to understand how to handle real truth from bwamem
    //{
    //    if(truth.is_open()) {
    //
    //        std::string read;
    //        int start;
    //        int end;
    //        while(truth >> read >> start >> end)
    //        {
    //            
    //            intervals.push_back(Interval<std::string>(start, end, read));
    //            queries.push_back(Interval<std::string>(start, end, read));
    //        }
    //
    //    } else std::cout << "Error opening the ground truth file" << endl;
    //}
    return trueOvls;
}

/* Compute the overlap length between potential overlapping reads pairs */
int computeLength(unordered_map<string,readInfo> & readMap, string & col_nametag, string & row_nametag) 
{

    int alignment_length = 0;

    unordered_map<string,readInfo>::const_iterator jit;
    unordered_map<string,readInfo>::const_iterator iit;

    jit = readMap.find(col_nametag); // col name
    iit = readMap.find(row_nametag); // row name 

    //cout << "ref 1: " << iit->second.ref << endl;
    //cout << "ref 2: " << jit->second.ref << endl;

    if(iit->second.ref == jit->second.ref)
    {
        //cout << "ref ==" << endl;
        if(iit->second.start < jit->second.start) {
            if(iit->second.end > jit->second.start)
            {
                alignment_length = min((iit->second.end - jit->second.start), (jit->second.end - jit->second.start));
            }
        }
        else if(iit->second.start > jit->second.start) 
        {
            if(jit->second.end > iit->second.start)
            {
                alignment_length = min((jit->second.end - iit->second.start), (iit->second.end - iit->second.start));
            }

        } 
        //else if(iit->second.start > jit->second.end) // to avoid R1.begin - R2. end >= THR and the length of R2 is >= THR
        //{
        //    alignment_length = 0;
        //}
        else
        {
            alignment_length = min((jit->second.end - iit->second.start), (iit->second.end - iit->second.start)); 
        }
    } //else cout << "!ref" << endl;
    //cout << alignment_length << endl;
    return alignment_length;
}

void benchmarkingAl(ifstream & groundtruth, ifstream & bella, ifstream & minimap, ifstream & mhap, 
    ifstream & blasr, ifstream & daligner, bool simulated, int minOvl, int nread) // add blasr && daligner && mhap && bella
{
    vector<vectInfo> truthInfo;
    unordered_map<string,readInfo> readMap;
    map<pair<string,string>, bool> checkBella;
    map<pair<string,string>, bool> checkMinimap;
    map<pair<string,string>, bool> checkMhap;
    map<pair<string,string>, bool> checkBlasr;
    map<pair<string,string>, bool> checkDaligner;
    map<pair<string,string>, bool>::iterator it;
    
    int alignment_length;
    
    double ovlsbella = 0, truebella = 0;
    double ovlsminimap = 0, trueminimap = 0;
    double ovlsmhap = 0, truemhap = 0;
    double ovlsblasr = 0, trueblasr = 0;
    double ovlsdal = 0, truedal = 0;

    cout << "\nbuilding the ground truth" << endl;

    if(simulated)
    {
        if(groundtruth.is_open())
        {
            string ref;
            string prev;
            string read;
            string dontcare1;
            string dontcare2;
            refInfo ovlInfo;
            readInfo perRead;
            vectInfo vectOf;
            int start;
            int end;

            while(groundtruth >> ref >> start >> end >> read >> dontcare1 >> dontcare2)
            {
                perRead.ref = ref;
                perRead.start = start;
                perRead.end = end;
                readMap.insert(make_pair("@"+read,perRead));

                //cout << "read: " << read << endl;

                ovlInfo.ref = ref;
                ovlInfo.read = "@" + read;
                ovlInfo.start = start;
                ovlInfo.end = end;
                vectOf.push_back(ovlInfo); // all the element of a chromosome

                if(ref != prev)
                {
                    truthInfo.push_back(vectOf); // its size should be the number of different references i.e. chromosomes
                    vectOf.clear();
                }

                prev = ref;
            }
            cout << "num reads: " << readMap.size() << endl;
            cout << "num chromosomes: " << truthInfo.size() << endl;
        } else cout << "Error opening the ground truth file" << endl;
    }
    //else // check how real data output is (sam file)
    //{
    //    if(groundtruth.is_open())
    //    {
    //        string read;
    //        int start;
    //        int end;
    //
    //        while(groundtruth >> read >> start >> end)
    //        {
    //            std::string nametag = "@" + read;
    //            int start_v = start;
    //            int end_v = end;
    //            coords = make_pair(start_v, end_v);
    //            ovlsmap.insert(std::pair<std::string, std::pair<int, int>>(nametag,coords));
    //        }
    //    } else std::cout << "Error opening the ground truth file" << endl;
    //}
    
    groundtruth.clear();
    groundtruth.seekg(0, ios::beg);

    cout << "computing BELLA recall/precision" << endl;
    if(bella.is_open())
    {
        string line;
        while(getline(bella, line))
        {
            ovlsbella++;
            //cout << "#" << ovlsbella << endl;
            stringstream lineStream(line);
            string col_nametag, row_nametag;

            getline(lineStream, col_nametag, ' ');
            getline(lineStream, row_nametag, ' ');
                                           // remove self aligned paired from bella output
            if(col_nametag == row_nametag) // to be sure to not count self aligned pairs
                ovlsbella--;
            else
            {
                it = checkBella.find(make_pair(col_nametag, row_nametag));
                if(it == checkBella.end())
                {       
                    checkBella.insert(make_pair(make_pair(col_nametag, row_nametag), true));
                    /* Compute the overlap length between potential overlapping reads pairs */
            
                    alignment_length = computeLength(readMap, col_nametag, row_nametag); // change function
                    if(alignment_length >= minOvl)
                        truebella++;
                        //cout << "aln: " << truebella << endl;
                }
            }
        }
    }

    cout << "computing Minimap recall/precision" << endl;
    if(minimap.is_open())
    {
        string line;
        while(getline(minimap, line))
        {
            ovlsminimap++;
            stringstream lineStream(line);
            string col_nametag, row_nametag, dontcare1, dontcare2, dontcare3, dontcare4;
    
            getline(lineStream, col_nametag, '\t' );
            getline(lineStream, dontcare1, '\t' );
            getline(lineStream, dontcare2, '\t' );
            getline(lineStream, dontcare3, '\t' );
            getline(lineStream, dontcare4, '\t' );
            getline(lineStream, row_nametag, '\t' );
    
            col_nametag = "@" + col_nametag;
            row_nametag = "@" + row_nametag;

            if(col_nametag == row_nametag) // to be sure to not count self aligned pairs
                ovlsminimap--;
            else
            {
                it = checkMinimap.find(make_pair(col_nametag, row_nametag));
                if(it == checkMinimap.end())
                {
                    checkMinimap.insert(make_pair(make_pair(col_nametag, row_nametag), true));
                    /* Compute the overlap length between potential overlapping reads pairs */
                    alignment_length = computeLength(readMap, col_nametag, row_nametag);
                    if(alignment_length >= minOvl)
                        trueminimap++;
                }
            }
        }
    }

    cout << "computing MHAP sensitive recall/precision" << endl;
    if(mhap.is_open())
    {
        string line;
        while(getline(mhap, line))
        {
            ovlsmhap++;
            stringstream lineStream(line);
            string col_nametag, row_nametag;

            getline(lineStream, col_nametag, ' ');
            getline(lineStream, row_nametag, ' ');

            col_nametag = "@" + col_nametag;
            row_nametag = "@" + row_nametag;

            if(col_nametag == row_nametag) // to be sure to not count self aligned pairs
                ovlsmhap--;
            else
            {
                it = checkMhap.find(make_pair(col_nametag, row_nametag));
                if(it == checkMhap.end())
                {
                    checkMhap.insert(make_pair(make_pair(col_nametag, row_nametag), true));
                    /* Compute the overlap length between potential overlapping reads pairs */
                    alignment_length = computeLength(readMap, col_nametag, row_nametag);
                    if(alignment_length >= minOvl)
                        truemhap++;
                }
            }
        }
    }

    cout << "computing BLASR recall/precision" << endl;
    if(blasr.is_open())
    {
        string line;
        while(getline(blasr, line))
        {
            ovlsblasr++;
            stringstream lineStream(line);
            string col_nametag, row_nametag;

            getline(lineStream, col_nametag, ' ');
            getline(lineStream, row_nametag, ' ');

            col_nametag = "@" + col_nametag;
            row_nametag = "@" + row_nametag;

            if(col_nametag == row_nametag) // to be sure to not count self aligned pairs
                ovlsblasr--;
            else
            {
                it = checkBlasr.find(make_pair(col_nametag, row_nametag));
                if(it == checkBlasr.end())
                {
                    checkBlasr.insert(make_pair(make_pair(col_nametag, row_nametag), true));
                    /* Compute the overlap length between potential overlapping reads pairs */
                    alignment_length = computeLength(readMap, col_nametag, row_nametag);
                    if(alignment_length >= minOvl)
                        trueblasr++;
                }
            }
        }
    }

    cout << "computing DALIGNER recall/precision" << endl;
    if(daligner.is_open())
    {
        string line;
        while(getline(daligner, line))
        {
            ovlsdal++;
            stringstream lineStream(line);
            string col_nametag, row_nametag;

            getline(lineStream, col_nametag, ' ');
            getline(lineStream, row_nametag, ' ');

            col_nametag = col_nametag;
            row_nametag = row_nametag;

            if(col_nametag == row_nametag) // to be sure to not count self aligned pairs
                ovlsblasr--;
            else
            {
                it = checkBlasr.find(make_pair(col_nametag, row_nametag));
                if(it == checkBlasr.end())
                {
                    checkBlasr.insert(make_pair(make_pair(col_nametag, row_nametag), true));
                    /* Compute the overlap length between potential overlapping reads pairs */
                    alignment_length = computeLength(readMap, col_nametag, row_nametag);
                    if(alignment_length >= minOvl)
                        trueblasr++;
                }
            }
        }
    }

    groundtruth.clear();
    groundtruth.seekg(0, ios::beg);
    double truetruth = trueOv(truthInfo, simulated, minOvl);

    // as -S option count overlaps only once 
    // (A ov B, but not B ov A), while all other 
    // (included the ground truth) count all ovls and the self-ovls
    double trueminimap_tuned = trueminimap*2;

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
