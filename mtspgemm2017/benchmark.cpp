#include "benchmark.h" 
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
#include <memory>

int main (int argc, char* argv[]) {

	// This program requires as input (so far) Minimap's overlaps + the ground truth computed with BWA-MEM
	// Later Daligner's, Blasr's and Mhap's overlaps will be required as well
	if(argc < 7)
    {
        cout << "not enough parameters, usage: " << endl;
        cout << "./benchmark <ground-truth> out.bella out.mini out.mhap out.blasr out.daligner" << endl;
    }

    std::ifstream ground(argv[1]);
	std::ifstream bella(argv[2]);
    std::ifstream minimap(argv[3]);
    std::ifstream mhap(argv[4]);
    // std::ifstream mhap_al(argv[5]);
    std::ifstream blasr(argv[5]);
    std::ifstream daligner(argv[6]);
    // std::ifstream daligner("out.m4");

    // getMetrics(bella, mhap, blasr);  
    benchmarkingAl(ground, bella, minimap, mhap, blasr, daligner);  

}