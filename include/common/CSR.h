#ifndef _CSR_H_
#define _CSR_H_

#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include "Deleter.h"
#include "CSC.h"
#include "Triple.h"

using namespace std;

template <class IT, class NT>
class CSR
{ 
public:
    CSR():nnz(0), rows(0), cols(0),zerobased(true) {}
	CSR(IT mynnz, IT m, IT n):nnz(mynnz),rows(m),cols(n),zerobased(true)
	{
        // Constructing empty Csc objects (size = 0) are allowed (why wouldn't they?).
        assert(rows != 0);
        
		rowptr = new IT[rows+1];
		if(nnz > 0)
		{
			colids = new IT[nnz];
			values = new NT[nnz];
		}
	}
    CSR (const CSC<IT,NT> & csc);   // CSC -> CSR conversion
    CSR (const CSR<IT,NT> & rhs);	// copy constructor
	CSR<IT,NT> & operator=(const CSR<IT,NT> & rhs);	// assignment operator
    
    ~CSR()
	{
        if( nnz > 0 )
            DeleteAll(colids, values);
        if( rows > 0 )
            delete [] rowptr;
	}
    bool ConvertOneBased()
    {
	if(!zerobased)	// already one-based
		return false; 
	transform(rowptr, rowptr + rows + 1, rowptr, bind2nd(plus<IT>(), static_cast<IT>(1)));
	transform(colids, colids + nnz, colids, bind2nd(plus<IT>(), static_cast<IT>(1)));
	zerobased = false;
	return true;
    }
    bool isEmpty()
    {
        return ( nnz == 0 );
    }
    void Sorted();
    
	IT rows;	
	IT cols;
	IT nnz; // number of nonzeros
    
    IT * rowptr;
    IT * colids;
    NT * values;
    bool zerobased;
};

#include "CSR.cpp"
#endif
