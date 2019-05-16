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

// overlap estimation from start/end alignment according to the overlapper's output
int estimeOv(int cstart, int cend, int clen, int rstart, int rend, int rlen)
{
    int diffc = cend-cstart;
    int diffr = rend-rstart;
    int left = min(cstart, rstart);
    int right = min(clen-cend, rlen-rend);

    int ov = left+right+(diffc+diffr)/2;
    return ov;
}

// from mecat index file to a map<uint32_t,string> : map<index,read-name>
void file2map(ifstream & idx2read, mecat_m & namestable)
{
    string num, name, seq;
    uint32_t idx;

    if(idx2read.is_open())
    {   
        string line;
        while(getline(idx2read, line))
        {
            stringstream lineStream(line);

            getline(lineStream, num, ' ');
            getline(lineStream, name, ' ' );

            getline(idx2read, seq); // sequence on new line

            idx = stoi(num);
            name.erase(0,1); // remove first char'>'
            namestable.insert(make_pair(idx, name));
        }
        cout << "MECAT idx2read table created." << endl;
    }
    else
    {
        cout << "Error creating names table from idx2read." << endl;
        exit(1);
    }
}

// from mecat numeric id to read name
string idx2read(uint32_t idx, mecat_m & namestable)
{
    mecat_m::iterator it;
    string name;

    it = namestable.find(idx);
    if(it != namestable.end())
    {
        return namestable[idx];
    }
    else
    {
        cout << "Read " << idx << " not present in MECAT output." << endl;
        exit(1);
    }
}


// mmap_m is a multimap where the read name is the key and the value is a 
// data struct containing reference and read name, and start and end alignment 
// of the read in the reference genome
// cname and rname are the (name of) pair of sequences from the overlapper's output
// this funtion checks if a pair of sequences is a true overlap according to the ground truth (reads mapped to reference)
int computeLength(mmap_m & seqmap, string & cname, string & rname) 
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