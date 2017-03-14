#ifndef _UTILITY_H
#define _UTILITY_H

#include <stdlib.h>
#include <stdint.h>
#include <climits>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;
#define 	EPSILON   0.001


template <class T>
struct ErrorTolerantEqual:
public binary_function< T, T, bool >
{
   ErrorTolerantEqual(const T & myepsilon):epsilon(myepsilon) {};
   inline bool operator() (const T & a, const T & b) const
   {
   	// According to the IEEE 754 standard, negative zero and positive zero should 
   	// compare as equal with the usual (numerical) comparison operators, like the == operators of C++ 
   
  	if(a == b)      // covers the "division by zero" case as well: max(a,b) can't be zero if it fails
   		return true;    // covered the integral numbers case
   
   	return ( std::abs(a - b) < epsilon || (std::abs(a - b) / max(std::abs(a), std::abs(b))) < epsilon ) ; 
   }
   T epsilon;
};

// Because identify reports ambiguity in PGI compilers
template<typename T>
struct myidentity : public std::unary_function<T, T>
{
    const T operator()(const T& x) const
    {
        return x;
    }
};

template<typename _ForwardIterator, typename _StrictWeakOrdering>
bool my_is_sorted(_ForwardIterator __first, _ForwardIterator __last,  _StrictWeakOrdering __comp)
{
   if (__first == __last)
   	return true;
   
   _ForwardIterator __next = __first;
   for (++__next; __next != __last; __first = __next, ++__next)
   	if (__comp(*__next, *__first))
   		return false;
   	return true;
};

template <typename ITYPE>
ITYPE CumulativeSum (ITYPE * arr, ITYPE size)
{
    ITYPE prev;
    ITYPE tempnz = 0 ;
    for (ITYPE i = 0 ; i < size ; ++i)
    {
		prev = arr[i];
		arr[i] = tempnz;
		tempnz += prev ;
    }
    return (tempnz) ;		    // return sum
}


template<typename _ForwardIter, typename T>
void iota(_ForwardIter __first, _ForwardIter __last, T __value)
{
	while (__first != __last)
     		*__first++ = __value++;
}
	
template<typename T, typename I>
T ** allocate2D(I m, I n)
{
	T ** array = new T*[m];
	for(I i = 0; i<m; ++i) 
		array[i] = new T[n];
	return array;
}

template<typename T, typename I>
void deallocate2D(T ** array, I m)
{
	for(I i = 0; i<m; ++i) 
		delete [] array[i];
	delete [] array;
}


template < typename T >
struct absdiff : binary_function<T, T, T>
{
        T operator () ( T const &arg1, T const &arg2 ) const
        {
                using std::abs;
                return abs( arg1 - arg2 );
        }
};

/* This function will return n % d.
   d must be one of: 1, 2, 4, 8, 16, 32, … */
inline unsigned int getModulo(unsigned int n, unsigned int d)
{
	return ( n & (d-1) );
} 

// Same requirement (d=2^k) here as well
inline unsigned int getDivident(unsigned int n, unsigned int d)
{
	while((d = d >> 1))
		n = n >> 1;
	return n;
}

#endif

