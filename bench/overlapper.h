#include "utils.h"
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

void metricsBella(ifstream & th, ifstream & bf, bool sim, int minOv, string & outfile, mmap_m & seqmap, uint32_t numth)
{
    bella_m metamap;
    bella_m::iterator it;
    int alignmentLength, ov;
    string cname, rname, nkmer, score, rev, cstart, cend, clen, rstart, rend, rlen, ovlen;
    uint32_t tpBella = 0, fpBella = 0, fnBella = 0, tnBella = 0;

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
            getline(lnstream, nkmer, '\t');
            getline(lnstream, score, '\t');
            getline(lnstream, ovlen, '\t');
            getline(lnstream, rev, '\t');
            getline(lnstream, cstart, '\t');
            getline(lnstream, cend, '\t');
            getline(lnstream, clen, '\t');
            getline(lnstream, rstart, '\t');
            getline(lnstream, rend, '\t');
            getline(lnstream, rlen, '\t');

            ov = stoi(ovlen);

            cname = "@" + cname;
            rname = "@" + rname;

            if(cname != rname) // not count self-pair
            {   // Check for BELLA outputting more than 1 alignment/pair
                it = metamap.find(make_pair(cname,rname)); // this shouldn't be useful
                if(it == metamap.end())
                {
                    if(ov >= minOv)
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength >= minOv) // TP
                        {
                            tpBella++;
    
                            metadata.st = 0;
                            metadata.cstart = cstart;
                            metadata.cend = cend;
                            metadata.clen = clen;
                            metadata.rstart = rstart;
                            metadata.rend = rend;
                            metadata.rlen = rlen;
                            metadata.nkmer = nkmer;
                            metadata.score = score;
                            metadata.ov = ov;
                            metadata.rev = rev;
                            metamap.insert(make_pair(make_pair(cname,rname),metadata));
    
                        }
                        else // FP
                        {
                            fpBella++;
    
                            metadata.st = 1;
                            metadata.cstart = cstart;
                            metadata.cend = cend;
                            metadata.clen = clen;
                            metadata.rstart = rstart;
                            metadata.rend = rend;
                            metadata.rlen = rlen;
                            metadata.nkmer = nkmer;
                            metadata.score = score;
                            metadata.ov = ov;
                            metadata.rev = rev;
                            metamap.insert(make_pair(make_pair(cname,rname),metadata));
    
                        }
                    }
                    else
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength >= minOv) // FP
                        {
                            fnBella++;
    
                            metadata.st = 1;
                            metadata.cstart = cstart;
                            metadata.cend = cend;
                            metadata.clen = clen;
                            metadata.rstart = rstart;
                            metadata.rend = rend;
                            metadata.rlen = rlen;
                            metadata.nkmer = nkmer;
                            metadata.score = score;
                            metadata.ov = ov;
                            metadata.rev = rev;
                            metamap.insert(make_pair(make_pair(cname,rname),metadata));
    
                        }
                        else // TN
                        {
                            tnBella++;
    
                            metadata.st = 2;
                            metadata.cstart = cstart;
                            metadata.cend = cend;
                            metadata.clen = clen;
                            metadata.rstart = rstart;
                            metadata.rend = rend;
                            metadata.rlen = rlen;
                            metadata.nkmer = nkmer;
                            metadata.score = score;
                            metadata.ov = ov;
                            metadata.rev = rev;
                            metamap.insert(make_pair(make_pair(cname,rname),metadata));
    
                        }
                    }
                }
                else if (it->second.st == 1) // tool claimed it as < 2kb the first time, but a multi mapping could result in a TP
                {
                    if(ov > minOv-1)
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpBella++;
                            fnBella--;

                            metadata.st = 0;
                            metadata.cstart = cstart;
                            metadata.cend = cend;
                            metadata.clen = clen;
                            metadata.rstart = rstart;
                            metadata.rend = rend;
                            metadata.rlen = rlen;
                            metadata.nkmer = nkmer;
                            metadata.score = score;
                            metadata.ov = ov;
                            metadata.rev = rev;
                            metamap[make_pair(cname,rname)] = metadata;
                        }
                    }
                }
                else if (it->second.st == 2) // tool claimed it as < 2kb the first time, but a multi mapping could result in a FP
                {
                    if(ov > minOv-1)
                    {
                        fpBella++;
                        tnBella--;

                        metadata.st = 1;
                        metadata.cstart = cstart;
                        metadata.cend = cend;
                        metadata.clen = clen;
                        metadata.rstart = rstart;
                        metadata.rend = rend;
                        metadata.rlen = rlen;
                        metadata.nkmer = nkmer;
                        metadata.score = score;
                        metadata.ov = ov;
                        metadata.rev = rev;
                        metamap[make_pair(cname,rname)] = metadata;
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

    ofstream output(outfile);
    if(output.is_open())
        for(it = metamap.begin(); it != metamap.end(); ++it)
        {
                output << it->second.st << ',' << it->second.score << ',' << it->second.ov << ',' << it->second.nkmer << ',' << it->first.first << ',' << it->first.second << ',' << it->second.cstart
                    << ',' << it->second.cend << ',' << it->second.clen << ',' << it->second.rstart << ',' << it->second.rend << ',' << it->second.rlen << endl;
        }
    output.close();

    cout << "BELLA:" << endl;
    cout << "   .Pairs = " << tpBella+fpBella+tnBella+fnBella << endl;
    cout << "   .FP = " << fpBella << endl;
    cout << "   .TP = " << tpBella << endl;
    cout << "   .TN = " << tnBella << endl;
    cout << "   .Recall = " << ((double)(tpBella*2)/(double)(numth+fnBella))*100 << "%" << endl;
    cout << "   .Precision = " << ((double)(tpBella)/(double)(tpBella+fpBella))*100 << "%" << endl;
    cout << "   .Specificity = " << ((double)(tnBella)/(double)(tnBella+fpBella))*100 << "%" << "\n" << endl;

    bf.close();
}

void metricsBellaPAF(ifstream & th, ifstream & bf, bool sim, int minOv, string & outfile, mmap_m & seqmap, uint32_t numth)
{
    bella_m metamap;
    bella_m::iterator it;
    int alignmentLength, ov;
    string cname, rname, score, rev, cstart, cend, clen, rstart, rend, rlen, ovlen, mapqv;
    uint32_t tpBella = 0, fpBella = 0, tnBella = 0;

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
            getline(lnstream, clen, '\t');
            getline(lnstream, cstart, '\t');
            getline(lnstream, cend, '\t');
            getline(lnstream, rev, '\t');
            getline(lnstream, rname, '\t');
            getline(lnstream, rlen, '\t');
            getline(lnstream, rstart, '\t');
            getline(lnstream, rend, '\t');
            getline(lnstream, score, '\t');
            getline(lnstream, ovlen, '\t');
            getline(lnstream, mapqv, '\t'); // missing and equal to 255     

            ov = stoi(ovlen);

            cname = "@" + cname;
            rname = "@" + rname;

            if(cname != rname) // not count self-pair
            {   // Check for BELLA outputting more than 1 alignment/pair
                it = metamap.find(make_pair(cname,rname)); // this shouldn't be useful
                if(it == metamap.end())
                {
                    if(ov >= minOv)
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength >= minOv) // TP
                        {
                            tpBella++;
    
                            metadata.st = 0;
                            metadata.cstart = cstart;
                            metadata.cend = cend;
                            metadata.clen = clen;
                            metadata.rstart = rstart;
                            metadata.rend = rend;
                            metadata.rlen = rlen;
                            metadata.score = score;
                            metadata.ov = ov;
                            metadata.rev = rev;
                            metamap.insert(make_pair(make_pair(cname,rname),metadata));
    
                        }
                        else // FP
                        {
                            fpBella++;
    
                            metadata.st = 1;
                            metadata.cstart = cstart;
                            metadata.cend = cend;
                            metadata.clen = clen;
                            metadata.rstart = rstart;
                            metadata.rend = rend;
                            metadata.rlen = rlen;
                            metadata.score = score;
                            metadata.ov = ov;
                            metadata.rev = rev;
                            metamap.insert(make_pair(make_pair(cname,rname),metadata));
    
                        }
                    }
                    else
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength >= minOv) // FP
                        {
                            fpBella++;
    
                            metadata.st = 1;
                            metadata.cstart = cstart;
                            metadata.cend = cend;
                            metadata.clen = clen;
                            metadata.rstart = rstart;
                            metadata.rend = rend;
                            metadata.rlen = rlen;
                            metadata.score = score;
                            metadata.ov = ov;
                            metadata.rev = rev;
                            metamap.insert(make_pair(make_pair(cname,rname),metadata));
    
                        }
                        else // TN
                        {
                            tnBella++;
    
                            metadata.st = 2;
                            metadata.cstart = cstart;
                            metadata.cend = cend;
                            metadata.clen = clen;
                            metadata.rstart = rstart;
                            metadata.rend = rend;
                            metadata.rlen = rlen;
                            metadata.score = score;
                            metadata.ov = ov;
                            metadata.rev = rev;
                            metamap.insert(make_pair(make_pair(cname,rname),metadata));
    
                        }
                    }
                }
                else if (it->second.st == 1) // tool claimed it as < 2kb the first time, but a multi mapping could result in a TP
                {
                    if(ov > minOv-1)
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpBella++;
                            fpBella--;

                            metadata.st = 0;
                            metadata.cstart = cstart;
                            metadata.cend = cend;
                            metadata.clen = clen;
                            metadata.rstart = rstart;
                            metadata.rend = rend;
                            metadata.rlen = rlen;
                            metadata.score = score;
                            metadata.ov = ov;
                            metadata.rev = rev;
                            metamap[make_pair(cname,rname)] = metadata;
                        }
                    }
                }
                else if (it->second.st == 2) // tool claimed it as < 2kb the first time, but a multi mapping could result in a FP
                {
                    if(ov > minOv-1)
                    {
                        fpBella++;
                        tnBella--;

                        metadata.st = 1;
                        metadata.cstart = cstart;
                        metadata.cend = cend;
                        metadata.clen = clen;
                        metadata.rstart = rstart;
                        metadata.rend = rend;
                        metadata.rlen = rlen;
                        metadata.score = score;
                        metadata.ov = ov;
                        metadata.rev = rev;
                        metamap[make_pair(cname,rname)] = metadata;
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

    ofstream output(outfile);
    if(output.is_open())
        for(it = metamap.begin(); it != metamap.end(); ++it)
        {
                output << it->second.st << ',' << it->second.score << ',' << it->second.ov << ',' << it->first.first << ',' << it->first.second << ',' << it->second.cstart
                    << ',' << it->second.cend << ',' << it->second.clen << ',' << it->second.rstart << ',' << it->second.rend << ',' << it->second.rlen << endl;
        }
    output.close();

    cout << "BELLA:" << endl;
    cout << "   .Pairs = " << tpBella+fpBella+tnBella << endl;
    cout << "   .FP = " << fpBella << endl;
    cout << "   .TP = " << tpBella << endl;
    cout << "   .TN = " << tnBella << endl;
    cout << "   .Recall = " << ((double)(tpBella*2)/(double)(numth))*100 << "%" << endl;
    cout << "   .Precision = " << ((double)(tpBella)/(double)(tpBella+fpBella))*100 << "%" << endl;
    cout << "   .Specificity = " << ((double)(tnBella)/(double)(tnBella+fpBella))*100 << "%" << "\n" << endl;

    bf.close();
}

void metricsMinimap2(ifstream & th, ifstream & mf, bool sim, int minOv, mmap_m & seqmap, uint32_t numth)
{
    check_m inMinimap;
    check_m::iterator it;

    int alignmentLength, ov;
    string cname, rname, nkmer, score, rev, cstart, cend, clen, rstart, rend, rlen;
    size_t tpMinimap = 0, fpMinimap = 0, tnMinimap = 0;

    /* st = 0: tp */
    /* st = 1: fp */
    /* st = 2: tn */

    cout << "MINIMAP2 evaluation started ..." << endl;
    if(mf.is_open())
    {
        string ln;
        while(getline(mf, ln))
        {
            stringstream lineStream(ln);

            getline(lineStream, cname, '\t' );
            getline(lineStream, clen, '\t' );
            getline(lineStream, cstart, '\t' );
            getline(lineStream, cend, '\t' );
            getline(lineStream, rev, '\t' );
            getline(lineStream, rname, '\t' );
            getline(lineStream, rlen, '\t' );
            getline(lineStream, rstart, '\t' );
            getline(lineStream, rend, '\t' );

            cname = "@" + cname;
            rname = "@" + rname;

            if(cname != rname) // not count self-pair
            {
                it = inMinimap.find(make_pair(cname,rname)); // this shouldn't be useful
                if(it == inMinimap.end())
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpMinimap++;
                            inMinimap.insert(make_pair(make_pair(cname, rname), 0));
                        }
                        else // FP
                        {
                            fpMinimap++;
                            inMinimap.insert(make_pair(make_pair(cname, rname), 1));
                        }
                    }
                    else
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // FP
                        {
                            fpMinimap++;
                            inMinimap.insert(make_pair(make_pair(cname, rname), 1));
                        }
                        else // tn
                        {
                            tnMinimap++;
                            inMinimap.insert(make_pair(make_pair(cname, rname), 2));
                        }
                    }
                }
                else if (it->second == 1) // tool claimed it as < 2kb the first time, but a multi mapping could result in a TP
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpMinimap++;
                            fpMinimap--;
                            inMinimap[make_pair(cname,rname)] = 0;
                        }
                    }
                }
                else if (it->second == 2) // tool claimed it as < 2kb the first time, but a multi mapping could result in a FP
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        fpMinimap++;
                        tnMinimap--;
                        inMinimap[make_pair(cname,rname)] = 1;
                    }
                }
            }
        }
    }
    else
    {
        cout << "Error opening MINIMAP2 output" << endl;
        exit(1);
    }

    // as -S option count overlaps only once (A ov B, but not B ov A), while all other (included the ground truth) count all ovls and the self-ovls
    cout << "MINIMAP2:" << endl;
    cout << "   .Pairs = " << tpMinimap+fpMinimap+tnMinimap << endl;
    cout << "   .FP = " << fpMinimap << endl;
    cout << "   .TP = " << tpMinimap << endl;
    cout << "   .TN = " << tnMinimap << endl;
    cout << "   .Recall = " << ((double)(tpMinimap*2)/(double)(numth))*100 << "%" << endl;
    cout << "   .Precision = " << ((double)(tpMinimap)/(double)(tpMinimap+fpMinimap))*100 << "%" << endl;
    cout << "   .Specificity = " << ((double)(tnMinimap)/(double)(tnMinimap+fpMinimap))*100 << "%" << "\n" << endl;

    mf.close();
}

void metricsMhap(ifstream & th, ifstream & pf, bool sim, int minOv, mmap_m & seqmap, uint32_t numth)
{
    check_m inMhap;
    check_m::iterator it;
    int alignmentLength, ov;
    string cname, rname, nkmer, score, rev, cstart, cend, clen, rstart, rend, rlen;
    uint32_t tpMhap = 0, fpMhap = 0, tnMhap = 0;

    /* st = 0: tp */
    /* st = 1: fp */
    /* st = 2: tn */

    cout << "MHAP evaluation started ..." << endl;
    if(pf.is_open())
    {
        string err, cstr, rstr;
        while(pf >> cname >> rname >> err >> nkmer >> cstr >> cstart >> cend >> clen >> rstr >> rstart >> rend >> rlen)
        {
            cname = "@" + cname;
            rname = "@" + rname;
            if(cname != rname) // not count self-pair
            {
                it = inMhap.find(make_pair(cname,rname)); // output could contain multiple entries for a given pair
                if(it == inMhap.end())
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpMhap++;
                            inMhap.insert(make_pair(make_pair(cname, rname), 0));
                        }
                        else // FP
                        {
                            fpMhap++;
                            inMhap.insert(make_pair(make_pair(cname, rname), 1));
                        }
                    }
                    else
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // FP
                        {
                            fpMhap++;
                            inMhap.insert(make_pair(make_pair(cname, rname), 1));
                        }
                        else // tn
                        {
                            tnMhap++;
                            inMhap.insert(make_pair(make_pair(cname, rname), 2));
                        }
                    }
                }
                else if (it->second == 1) // tool claimed it as < 2kb the first time, but a multi mapping could result in a TP
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpMhap++;
                            fpMhap--;
                            inMhap[make_pair(cname,rname)] = 0;
                        }
                    }
                }
                else if (it->second == 2) // tool claimed it as < 2kb the first time, but a multi mapping could result in a FP
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        fpMhap++;
                        tnMhap--;
                        inMhap[make_pair(cname,rname)] = 1;
                    }
                }
            }
        }
    }
    else
    {
        cout << "Error opening MHAP output" << endl;
        exit(1);
    }

    cout << "MHAP:" << endl;
    cout << "   .Pairs = " << tpMhap+fpMhap+tnMhap << endl;
    cout << "   .FP = " << fpMhap << endl;
    cout << "   .TP = " << tpMhap << endl;
    cout << "   .TN = " << tnMhap << endl;
    cout << "   .Recall = " << ((double)(tpMhap)/(double)(numth))*100 << "%" << endl;
    cout << "   .Precision = " << ((double)(tpMhap)/(double)(tpMhap+fpMhap))*100 << "%" << endl;
    cout << "   .Specificity = " << ((double)(tnMhap)/(double)(tnMhap+fpMhap))*100 << "%" << "\n" << endl;

    pf.close();
}

void metricsMecat(ifstream & th, ifstream & cf, bool sim, int minOv, mmap_m & seqmap, uint32_t numth, ifstream & readnames)
{
    check_m inMecat;
    check_m::iterator it;
    int alignmentLength, ov;
    uint32_t tpMecat = 0, fpMecat = 0, tnMecat = 0;
    mecat_m namestable;

    file2map(readnames,namestable);

    /* st = 0: tp */
    /* st = 1: fp */
    /* st = 2: tn */

    cout << "MECAT evaluation started ..." << endl;
    if(cf.is_open())
    {
        string ln;
        while(getline(cf, ln))
        {
            stringstream lineStream(ln);
            string cid, rid, cname, rname, score, rev, cstart, cend, clen, rstart, rend, rlen, id, cstr, rstr;

            getline(lineStream, cid, '\t' );
            getline(lineStream, rid, '\t' );
            getline(lineStream, score, '\t' );
            getline(lineStream, id, '\t' );
            getline(lineStream, cstr, '\t' );
            getline(lineStream, cstart, '\t' );
            getline(lineStream, cend, '\t' );
            getline(lineStream, clen, '\t' );
            getline(lineStream, rstr, '\t' );
            getline(lineStream, rstart, '\t' );
            getline(lineStream, rend, '\t' );
            getline(lineStream, rlen, '\t' );

            cname = idx2read(stoi(cid), namestable);
            rname = idx2read(stoi(rid), namestable);

            cname = "@" + cname;
            rname = "@" + rname;

            if(cname != rname) // not count self-pair
            {
                it = inMecat.find(make_pair(cname,rname)); // output could contain multiple entries for a given pair
                if(it == inMecat.end())
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        alignmentLength = computeLength(seqmap,cname,rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpMecat++;
                            inMecat.insert(make_pair(make_pair(cname,rname), 0));
                        }
                        else // FP
                        {
                            fpMecat++;
                            inMecat.insert(make_pair(make_pair(cname,rname), 1));
                        }
                    }
                    else
                    {
                        alignmentLength = computeLength(seqmap,cname,rname);
                        if(alignmentLength > minOv-1) // FP
                        {
                            fpMecat++;
                            inMecat.insert(make_pair(make_pair(cname,rname), 1));
                        }
                        else // tn
                        {
                            tnMecat++;
                            inMecat.insert(make_pair(make_pair(cname,rname), 2));
                        }
                    }
                }
                else if (it->second == 1) // tool claimed it as < 2kb the first time, but a multi mapping could result in a TP
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        alignmentLength = computeLength(seqmap,cname,rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpMecat++;
                            fpMecat--;
                            inMecat[make_pair(cname,rname)] = 0;
                        }
                    }
                }
                else if (it->second == 2) // tool claimed it as < 2kb the first time, but a multi mapping could result in a FP
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        fpMecat++;
                        tnMecat--;
                        inMecat[make_pair(cname,rname)] = 1;
                    }
                }
            }
        }
    }
    else
    {
        cout << "Error opening MECAT output" << endl;
        exit(1);
    }

    cout << "MECAT:" << endl;
    cout << "   .Pairs = " << tpMecat+fpMecat+tnMecat << endl;
    cout << "   .FP = " << fpMecat << endl;
    cout << "   .TP = " << tpMecat << endl;
    cout << "   .TN = " << tnMecat << endl;
    cout << "   .Recall = " << ((double)(tpMecat*2)/(double)(numth))*100 << "%" << endl;
    cout << "   .Precision = " << ((double)(tpMecat)/(double)(tpMecat+fpMecat))*100 << "%" << endl;
    cout << "   .Specificity = " << ((double)(tnMecat)/(double)(tnMecat+fpMecat))*100 << "%" << "\n" << endl;

    cf.close();
}

void metricsBlasr(ifstream & th, ifstream & lf, bool sim, int minOv, mmap_m & seqmap, uint32_t numth)
{

    check_m inBlasr;
    check_m::iterator it;
    int alignmentLength, ov;
    string cname, rname, score, rev, cstart, cend, clen, rstart, rend, rlen, id, cstr, rstr, qv;
    uint32_t tpBlasr = 0, fpBlasr = 0, tnBlasr = 0;

    /* st = 0: tp */
    /* st = 1: fp */
    /* st = 2: tn */

    cout << "BLASR evaluation started ..." << endl;
    if(lf.is_open())
    {
        while(lf >> cname >> rname >> score >> id >> cstr >> cstart >> cend >> clen >> rstr >> rstart >> rend >> rlen >> qv)
        {
            cname = "@" + cname;
            rname = "@" + rname;
            if(cname != rname) // not count self-pair
            {
                it = inBlasr.find(make_pair(cname,rname)); // output could contain multiple entries for a given pair
                if(it == inBlasr.end())
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpBlasr++;
                            inBlasr.insert(make_pair(make_pair(cname, rname), 0));
                        }
                        else // FP
                        {
                            fpBlasr++;
                            inBlasr.insert(make_pair(make_pair(cname, rname), 1));
                        }
                    }
                    else
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // FP
                        {
                            fpBlasr++;
                            inBlasr.insert(make_pair(make_pair(cname, rname), 1));
                        }
                        else // tn
                        {
                            tnBlasr++;
                            inBlasr.insert(make_pair(make_pair(cname, rname), 2));
                        }
                    }
                }
                else if (it->second == 1) // tool claimed it as < 2kb the first time, but a multi mapping could result in a TP
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpBlasr++;
                            fpBlasr--;
                            inBlasr[make_pair(cname,rname)] = 0;
                        }
                    }
                }
                else if (it->second == 2) // tool claimed it as < 2kb the first time, but a multi mapping could result in a FP
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), stoi(rstart), stoi(rend), stoi(rlen));
                    if(ov > minOv-1)
                    {
                        fpBlasr++;
                        tnBlasr--;
                        inBlasr[make_pair(cname,rname)] = 1;
                    }
                }
            }
        }
    }
    else
    {
        cout << "Error opening BLASR output" << endl;
        exit(1);
    }

    cout << "BLASR:" << endl;
    cout << "   .Pairs = " << tpBlasr+fpBlasr+tnBlasr << endl;
    cout << "   .FP = " << fpBlasr << endl;
    cout << "   .TP = " << tpBlasr << endl;
    cout << "   .TN = " << tnBlasr << endl;
    cout << "   .Recall = " << ((double)(tpBlasr)/(double)(numth))*100 << "%" << endl;
    cout << "   .Precision = " << ((double)(tpBlasr)/(double)(tpBlasr+fpBlasr))*100 << "%" << endl;
    cout << "   .Specificity = " << ((double)(tnBlasr)/(double)(tnBlasr+fpBlasr))*100 << "%" << "\n" << endl;

    lf.close();
}

void metricsDaligner(ifstream & th, ifstream & df, bool sim, int minOv, multimap<string,readInfo> & seqmap, size_t numth)
{
    check_m inDaligner;
    check_m::iterator it;
    int alignmentLength, ov;
    string cname, rname, nkmer, score, rev, cstart, cend, clen, rstart, rend, rlen;
    uint32_t tpDaligner = 0, fpDaligner = 0, tnDaligner = 0;

    /* st = 0: tp */
    /* st = 1: fp */
    /* st = 2: tn */

    cout << "DALIGNER evaluation started ..." << endl;
    if(df.is_open())
        while(df >> cname >> rname >> rev >> cstart >> cend >> clen >> rstart >> rend >> rlen)
        {
            uint32_t intrstart = stoi(rstart);
            uint32_t intrend = stoi(rend);
            uint32_t intrlen = stoi(rlen);

            if(rev.compare("c") == 0)
            {
                intrstart = intrlen-intrend-1;
                intrend = intrlen-intrstart-1;
            }

            cname = "@" + cname;
            rname = "@" + rname;

            if(cname != rname) // not count self-pair
            {
                it = inDaligner.find(make_pair(cname,rname));  // output could contain multiple entries for a given pair
                if(it == inDaligner.end())
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), intrstart, intrend, intrlen);
                    if(ov > minOv-1)
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpDaligner++;
                            inDaligner.insert(make_pair(make_pair(cname, rname), 0));
                        }
                        else // FP
                        {
                            fpDaligner++;
                            inDaligner.insert(make_pair(make_pair(cname, rname), 1));
                        }
                    }
                    else
                    {
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // FP
                        {
                            fpDaligner++;
                            inDaligner.insert(make_pair(make_pair(cname, rname), 1));
                        }
                        else // tn
                        {
                            tnDaligner++;
                            inDaligner.insert(make_pair(make_pair(cname, rname), 2));
                        }
                    }
                }
                else if (it->second == 1) // tool claimed it as < 2kb the first time, but a multi mapping could result in a TP
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), intrstart, intrend, intrlen);
                    if(ov > minOv-1)
                        alignmentLength = computeLength(seqmap, cname, rname);
                        if(alignmentLength > minOv-1) // TP
                        {
                            tpDaligner++;
                            fpDaligner--;
                            inDaligner[make_pair(cname,rname)] = 0;
                        }
                }
                else if (it->second == 2) // tool claimed it as < 2kb the first time, but a multi mapping could result in a FP
                {
                    ov = estimeOv(stoi(cstart), stoi(cend), stoi(clen), intrstart, intrend, intrlen);
                    if(ov > minOv-1)
                    {
                        fpDaligner++;
                        tnDaligner--;
                        inDaligner[make_pair(cname,rname)] = 1;
                    }
                }
            }
        }
    else
    {
        cout << "Error opening DALIGNER output" << endl;
        exit(1);
    }

    cout << "DALIGNER:" << endl;
    cout << "   .Pairs = " << tpDaligner+fpDaligner+tnDaligner << endl;
    cout << "   .FP = " << fpDaligner << endl;
    cout << "   .TP = " << tpDaligner << endl;
    cout << "   .TN = " << tnDaligner << endl;
    cout << "   .Recall = " << ((double)(tpDaligner)/(double)(numth))*100 << "%" << endl;
    cout << "   .Precision = " << ((double)(tpDaligner)/(double)(tpDaligner+fpDaligner))*100 << "%" << endl;
    cout << "   .Specificity = " << ((double)(tnDaligner)/(double)(tnDaligner+fpDaligner))*100 << "%" << "\n" << endl;

    df.close();
}
