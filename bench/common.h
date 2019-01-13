#ifndef _COMMON_H_
#define _COMMON_H_

#include <string>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <map>
#include "IntervalTree.h"

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

struct myStruct {
    int st;
    string cstart;
    string cend;
    string clen;
    string rstart;
    string rend;
    string rlen;
    string nkmer;
    string score;
    int ov;
    string rev;
};

typedef vector<refInfo> info_v;
typedef map<uint32_t,string> mecat_m;
typedef map<pair<string,string>,int> check_m;
typedef map<pair<string,string>,myStruct> bella_m;
typedef vector<Interval<std::string>> intervals_v;
typedef map<string,info_v> truth_m;
typedef multimap<string,readInfo> mmap_m;

#endif