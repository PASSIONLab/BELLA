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
// Factorial function
//

double factorial(double n)
{
    if(n > 1)
        return n * factorial(n - 1);
    else
        return 1;
}

//
// Upper bound selection function, the lower bound is fixed and equal to 2 so far
//
int rbounds(int d, double e, int k)
{
    double a,b,c;
    double probability = 1;
    int m = 2;

    while(probability >= MIN_PROB)
    {
        m++;
        a = factorial(d)/(factorial(m)*factorial(d-m)); // it's fine 
        b = pow(1-e,(m*k));
        c = pow(1-pow(1-e,k),(d-m));
        
        probability = a*b*c;
    }
    return (m-1);
}

