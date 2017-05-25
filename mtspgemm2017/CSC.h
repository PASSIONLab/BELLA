#ifndef _CSC_H_
#define _CSC_H_

#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <tuple>
#include "Deleter.h"
#include "HeapEntry.h"

#include "Triple.h"
extern "C" {
#include "GTgraph/R-MAT/graph.h"
}

using namespace std;

template <class IT, class NT> // IT, NT li dichiaro runtime (polimorfismo parametrico)s
class CSC
{
public:
    CSC():nnz(0), rows(0), cols(0) {}
    CSC(IT mynnz, IT m, IT n, int nt):nnz(mynnz),rows(m),cols(n) // costruttore di default 
    {
        // Constructing empty Csc objects (size = 0) are not allowed.
        assert(nnz != 0 && cols != 0);
        
        colptr = new IT[cols+1]; 
        rowids = new IT[nnz];
        values = new NT[nnz];
    }
    CSC (Triple<IT,NT> * triples, IT mynnz, IT m, IT n); // altro costruttore di default
    
    template <typename AddOperation>
    CSC (vector<tuple<IT,IT,NT>> & tuple, IT m, IT n, AddOperation addop); // costruttore
    
    template <typename AddOperation>
    void MergeDuplicates (AddOperation addop); // method

    CSC(graph & G);
    CSC (IT * ri, IT * ci, NT * val, IT mynnz, IT m, IT n);
    CSC (const CSC<IT,NT> & rhs);		// copy constructor
    CSC<IT,NT> & operator=(const CSC<IT,NT> & rhs);	// assignment operator
    bool operator==(const CSC<IT,NT> & rhs); // ridefinizione ==
    //template <typename FT> CSC<IT,FT> & operator=(const CSC<IT,FT> & rhs);	// assignment operator
    //template <typename FT> bool operator==(const CSC<IT,FT> & rhs); // ridefinizione ==
    
    ~CSC() // distruttore
    {
        if( nnz > 0 )
            DeleteAll(rowids, values);
        if( cols > 0 )
            delete [] colptr;
    }
  
    bool isEmpty()
    {
        return ( nnz == 0 );
    }
    void Sorted();
    CSC<IT,NT> SpRef (const vector<IT> & ri, const vector<IT> & ci);
    CSC<IT,NT> SpRef1 (const vector<IT> & ri, const vector<IT> & ci);
    CSC<IT,NT> SpRef2 (const IT* ri, const IT rilen, const IT* ci, const IT cilen);
    void intersect (const IT* rowids_in, const NT* values_in, const IT len_in,
                                const IT* ri, const IT len_ri,
                    IT* rowids_out, NT* values_out, IT* len_out);

    IT rows;
    IT cols;
    IT nnz; // number of nonzeros
    IT totalcols;   // for the parallel case
    
    IT * colptr;
    IT * rowids;
    NT * values;
};

#include "CSC.cpp"
#endif
