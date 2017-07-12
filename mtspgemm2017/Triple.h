#ifndef _TRIPLE_H_
#define _TRIPLE_H_

#include <iostream>
#include <functional>
using namespace std;

template <class IT, class NT>
struct Triple
{
	IT row;	//  row index
	IT col;	//  col index
	NT val;		//  value
	Triple(IT myrow, IT mycol, NT myval):row(myrow),col(mycol),val(myval) {};
	Triple():row(0),col(0),val(0) {};

};

#endif
