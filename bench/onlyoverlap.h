#include "common.h"
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
#include <typeinfo>

using namespace std;

// mmap_m is a multimap where the read name is the key and the value is a 
// data struct containing reference and read name, and start and end alignment 
// of the read in the reference genome
// cname and rname are the (name of) pair of sequences from the overlapper's output
// this funtion checks if a pair of sequences is a true overlap according to the ground truth (reads mapped to reference)
int computeLength2(mmap_m & seqmap, string & cname, string & rname) 
{
    int len = 0, max = 0;

    mmap_m::iterator ccheck = seqmap.find(cname); // col name
    mmap_m::iterator rcheck = seqmap.find(rname); // row name 

    if(ccheck != seqmap.end() && rcheck != seqmap.end()) // needed as handling real dataset the aligned reads in sam file could be != the original number of reads
    {
        // this returns a pair representing the range of elements with key equal to the one we want
        pair<mmap_m::iterator,mmap_m::iterator> crange = seqmap.equal_range(cname);
        pair<mmap_m::iterator,mmap_m::iterator> rrange = seqmap.equal_range(rname);
        // this terates over the range
        for (mmap_m::iterator itC = crange.first; itC != crange.second; itC++)
        {
            for (mmap_m::iterator itR = rrange.first; itR != rrange.second; itR++)
            {
                if(itC->second.ref == itR->second.ref)
                {   
                    if(itC->second.start < itR->second.start) 
                    {
                        if(itC->second.end > itR->second.start)
                            len = min((itC->second.end - itR->second.start), (itR->second.end - itR->second.start));
                    }
                else if(itC->second.start > itR->second.start) 
                {
                    if(itR->second.end > itC->second.start)
                        len = min((itR->second.end - itC->second.start), (itC->second.end - itC->second.start));
                } 
                else len = min((itR->second.end - itC->second.start), (itC->second.end - itC->second.start)); 
                }

                if(len > max)
                    max = len;
            }
        }
    }
    return max;
}

void myBella(ifstream & th, ifstream & bf, bool sim, int minOv, mmap_m & seqmap, uint32_t numth)
{
    check_m metamap;
    check_m::iterator it;
    string cname, rname;
    uint32_t tpBella = 0, fpBella = 0, tnBella = 0;
    int alignmentLength;

    /* st = 0: tp */
    /* st = 1: fp */
    /* st = 2: tn */

    cout << "BELLA evaluation started ..." << endl;
    if(bf.is_open())
    {
        string ln;
        while(getline(bf, ln))
        {
            stringstream lnstream(ln);
            myStruct metadata;

            getline(lnstream, cname, '\t');
            getline(lnstream, rname, '\t');

            if(cname != rname) // not count self-pair
            {   // Check for BELLA outputting more than 1 alignment/pair
                it = metamap.find(make_pair(cname,rname)); // this shouldn't be useful
                if(it == metamap.end())
                {
                    alignmentLength = computeLength2(seqmap, cname, rname);
                    if(alignmentLength >= minOv) // TP
                    {
                        tpBella++;
                        metamap.insert(make_pair(make_pair(cname,rname),0));
                    }
                    else // FP
                    {
                        fpBella++;
                        metamap.insert(make_pair(make_pair(cname,rname),1));
                    }
                }
            }// if(cname != rname)
        } // while(getline(bf, ln))
    }
    else
    {
        cout << "Error opening BELLA output" << endl;
        exit(1);
    }

    cout << "BELLA:" << endl;
    cout << "   .Pairs = " << tpBella+fpBella << endl;
    cout << "   .FP = " << fpBella << endl;
    cout << "   .TP = " << tpBella << endl;
    cout << "   .Recall = " << ((double)(tpBella*2)/(double)(numth))*100 << "%" << endl;
    cout << "   .Precision = " << ((double)(tpBella)/(double)(tpBella+fpBella))*100 << "%" << endl;

    bf.close();
}

void myDiBella(ifstream & th, ifstream & bf, bool sim, int minOv, mmap_m & seqmap, uint32_t numth)
{
    check_m metamap;
    check_m::iterator it;
    string cid, rid, cname, rname;
    uint32_t tpDiBella = 0, fpDiBella = 0, tnDiBella = 0;
    int alignmentLength;

    /* st = 0: tp */
    /* st = 1: fp */
    /* st = 2: tn */

    cout << "diBELLA evaluation started ..." << endl;
    if(bf.is_open())
    {
        string ln;
        while(getline(bf, ln))
        {
            stringstream lnstream(ln);
            myStruct metadata;

            getline(lnstream, cid, ' ');
            getline(lnstream, rid, ' ');
            getline(lnstream, cname, ' ');
            getline(lnstream, rname, ' ');

            if(cname != rname) // not count self-pair
            {   // Check for BELLA outputting more than 1 alignment/pair
                it = metamap.find(make_pair(cname,rname)); // this shouldn't be useful
                if(it == metamap.end())
                {
                    alignmentLength = computeLength2(seqmap, cname, rname);
                    if(alignmentLength >= minOv) // TP
                    {
                        tpDiBella++;
                        metamap.insert(make_pair(make_pair(cname,rname),0));
                    }
                    else // FP
                    {
                        fpDiBella++;
                        metamap.insert(make_pair(make_pair(cname,rname),1));
                    }
                }
            }// if(cname != rname)
        } // while(getline(bf, ln))
    }
    else
    {
        cout << "Error opening diBELLA output" << endl;
        exit(1);
    }

    cout << "diBELLA:" << endl;
    cout << "   .Pairs = " << tpDiBella+fpDiBella+tnDiBella << endl;
    cout << "   .FP = " << fpDiBella << endl;
    cout << "   .TP = " << tpDiBella << endl;
    cout << "   .Recall = " << ((double)(tpDiBella*2)/(double)(numth))*100 << "%" << endl;
    cout << "   .Precision = " << ((double)(tpDiBella)/(double)(tpDiBella+fpDiBella))*100 << "%" << endl;

    bf.close();
}