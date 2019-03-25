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
#include "bound.hpp"

/**
 * @brief factorial
 * @param n
 * @return
 */
long double factorial(double n)
{
    if(n > 1)
        return n * factorial(n - 1);
    else
        return 1;
}

/**
 * @brief rbounds does upper bound selection,
 * the lower bound is fixed and equal to 2
 * @param d is the coverage
 * @param e is the error rate
 * @param k is the k-mer length
 * @return upper bound
 */
int computeUpper(int d, double e, int k)
{
    long double a,b,c;
    double probability = 1;
    double cumsum = 0, prev;
    int m = d;

    while(cumsum < MIN_PROB)
    {
        a = factorial(d)/(factorial(m)*factorial(d-m)); // it's fine 
        b = pow(1-e,(m*k));
        c = pow(1-pow(1-e,k),(d-m));
        
        probability = a*b*c;
        cumsum = cumsum + probability;

        if(cumsum == prev && cumsum < MIN_PROB)
            break;
        --m;
        prev = cumsum;
    }
    return (m+1);
}

int computeLower(int d, double e, int k)
{
    long double a,b,c;
    double probability = 1;
    double cumsum = 0, prev;
    int m = 2;

    while(cumsum < MIN_PROB)
    {
        a = factorial(d)/(factorial(m)*factorial(d-m)); // it's fine 
        b = pow(1-e,(m*k));
        c = pow(1-pow(1-e,k),(d-m));
        
        probability = a*b*c;
        cumsum = cumsum + probability;

        if(cumsum == prev && cumsum < MIN_PROB)
            break;
        ++m;
        prev = cumsum;
    }

    if (m-1 < 2)
        return 2;
    else 
        return (m-1);
}




