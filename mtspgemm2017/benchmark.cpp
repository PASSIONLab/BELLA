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
	if(argc < 3)
    {
        cout << "not enough parameters, usage: " << endl;
        cout << "./benchmark groundtruth.txt out.mini" << endl;
    }

    std::ifstream ground(argv[1]);
	// std::ifstream bella(argv[2]);
    std::ifstream minimap(argv[2]);
    // std::ifstream blasr("out.m4");
    // std::ifstream mhap("out.m4");
    // std::ifstream daligner("out.m4");

    // getMetrics(bella, mhap, blasr);  
    benchmarkingAl(ground, minimap);  

}