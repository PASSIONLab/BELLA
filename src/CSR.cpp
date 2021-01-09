#include "CSR.h"
#include "utility.h"

// copy constructor
template <class IT, class NT>
CSR<IT,NT>::CSR (const CSR<IT,NT> & rhs): nnz(rhs.nnz), rows(rhs.rows), cols(rhs.cols),zerobased(rhs.zerobased)
{
	if(nnz > 0)
	{
		values = new NT[nnz];
		colids = new IT[nnz];
        copy(rhs.values, rhs.values+nnz, values);
        copy(rhs.colids, rhs.colids+nnz, colids);
	}
	if ( rows > 0)
	{
		rowptr	= new IT[rows+1];
        copy(rhs.rowptr, rhs.rowptr+rows+1, rowptr);
	}
}

template <class IT, class NT>
CSR<IT,NT> & CSR<IT,NT>::operator= (const CSR<IT,NT> & rhs)
{
	if(this != &rhs)		
	{
		if(nnz > 0)	// if the existing object is not empty
		{
			delete [] colids;   // empty it
			delete [] values;
		}
		if(rows > 0)
		{
			delete [] rowptr;
		}

		nnz	= rhs.nnz;
		rows = rhs.rows;
		cols = rhs.cols;
		zerobased = rhs.zerobased;
		if(rhs.nnz > 0)	// if the copied object is not empty
		{
			values = new NT[nnz];
			colids = new IT[nnz];
            copy(rhs.values, rhs.values+nnz, values);
            copy(rhs.colids, rhs.colids+nnz, colids);
		}
		if(rhs.cols > 0)
		{
			rowptr	= new IT[rows+1];
            copy(rhs.rowptr, rhs.rowptr+rows+1, rowptr);
		}
	}
	return *this;
}

//! Construct a CSR object from a CSC
//! Accepts only zero based CSC inputs
template <class IT, class NT>
CSR<IT,NT>::CSR(const CSC<IT,NT> & csc):nnz(csc.nnz), rows(csc.rows), cols(csc.cols),zerobased(true)
{
	rowptr = new IT[rows+1];
    	colids = new IT[nnz];
    	values = new NT[nnz];

	IT * work = new IT[rows];	// workspace
    	std::fill(work, work+rows, (IT) 0); // initilized to zero
   
    	for (IT k = 0 ; k < nnz ; ++k)
    	{
        	IT tmp =  csc.rowids[k];
        	work [ tmp ]++ ;		// row counts (i.e, w holds the "row difference array")
	}

	if(nnz > 0)
	{
		rowptr[rows] = CumulativeSum (work, rows);		// cumulative sum of w
       	 	copy(work, work+rows, rowptr);

		IT last;
        	for (IT i = 0; i < cols; ++i) 
        	{
       	     		for (IT j = csc.colptr[i]; j < csc.colptr[i+1] ; ++j)
            		{
				colids[ last = work[ csc.rowids[j] ]++ ]  = i ;
				values[last] = csc.values[j] ;
            		}
        	}
	}
	delete [] work;
}


// check if sorted within rows?
template <class IT, class NT>
void CSR<IT,NT>::Sorted()
{
	bool sorted = true;
	for(IT i=0; i< rows; ++i)
	{
		sorted &= my_is_sorted (colids + rowptr[i], colids + rowptr[i+1], std::less<IT>());
    }
	cout << "Sorted ? " << sorted << endl;
}
