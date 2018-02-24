#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <array>
#include <tuple>
#include <queue>
#include <memory>
#include <stack>
#include <functional>
#include <cstring>
#include <string.h>
#include <cassert>
#include <ios>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <math.h>
#include <map>
#include <omp.h>
#include "rbounds.hpp"

//
// Binomial coefficient function
//
int bincoef(int n, int k)
{
    int num = 1;
 
    // C(n, k) = C(n, n-k)
    if (k > n-k)
        k = n - k;
 
    // Calculate value of [n*(n-1) * --- * (n-k+1)]/[k*(k-1)* ---- *1]
    for (int i = 0; i < k; ++i)
    {
        num *= (n-i);
        num /= (i+1);
    }
    return num;
}

//
// Upper bound selection function, the lower bound is fixed and equal to 2 so far
//
int rbounds(int d, float e, int k)
{
	int m = 2; 			// multiplicity in the fastq(s)
	float uprob = 1.0; 	// max probability 

	while(uprob >= MIN_PROB) 
	{
		++m;
		uprob = bincoef(d,m)*pow(1-e,m*k)*pow(1-pow(1-e,k),d-m);
	}
	--m; 
	return m;
}

