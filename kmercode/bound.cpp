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

//	GG: parallel factorial
long double factorial(long double number)
{
	long double fac = 1;
#pragma omp parallel for reduction(*:fac)
	for(int n = 2; n <= (int)number; ++n)
		fac *= n;
	return fac;
}

//	GG: depth, error rate, and k-mer length
int computeUpper(int d, double e, int k)
{
	long double a, b, c;
	long double dfact = factorial(d);
	long double bbase = (1-e);
	long double cbase = (1-pow(bbase, k));
	long double probability = 1;
	long double sum = 0, prev;
	int m = d;

	while(sum < MINPROB)
	{
		a = dfact / (factorial(m) * factorial(d-m));
		b = pow(bbase, (m*k));
		c = pow(cbase, (d-m));

		probability = a * b * c;
		sum = sum + probability;

		if(sum == prev && sum < MINPROB)
			break;
		--m;
		prev = sum;
	}
	return (m+1);
}

//	GG: depth, error rate, and k-mer length
int computeLower(int d, double e, int k)
{
	long double a, b, c;
	long double dfact = factorial(d);
	long double bbase = (1-e);
	long double cbase = (1-pow(bbase, k));
	long double probability = 1;
	long double sum = 0, prev;
	int mymin = 2;
	int m = mymin;

	while(sum < MINPROB)
	{
		a = dfact / (factorial(m) * factorial(d-m));
		b = pow(bbase, (m*k));
		c = pow(cbase, (d-m));

		probability = a * b * c;
		sum = sum + probability;

		if(sum == prev && sum < MINPROB)
			break;
		++m;
		prev = sum;
	}

	return std::max(m-1, mymin);
}




