#ifndef _TRANSPOSE_H_
#define _TRANSPOSE_H_
#include <atomic>

/** 
 ** while the name suggests CSR -> CSC, 
 ** one can simply use this to transpose a matrix in CSR format to its transpose also in CSR format
 ** or vice versa: transpose a matrix in CSC format to its transpose also in CSC format
 ** just read CSR related params as the "old matrix" and CSC related params as the "output/new matrix"
 ** currently only the nonsorted output case is implemented (i.e. rowids in each column of the output matrix is not sorted)
 **/
template <class IT, class NT>
void csr2csc_atomic_nosort (IT m, IT n, IT nnz, IT * csrRowPtr, IT * csrColIdx, NT * csrVal,
			IT * cscColPtr, IT * cscRowIdx, NT * cscVal)
{
	std::atomic<int> * atomicColPtr = new std::atomic<int>[n];
	for (IT i=0; i < n; i++)
	{
		atomicColPtr[i] = 0;
	}

	// construct an array of size nnz to record the relative
	// position of a nonzero element in corresponding column
	// this is what allows us to parallelize the last loop
	IT * dloc = new IT[nnz]();
	#pragma omp parallel for schedule(dynamic)
	for (IT i=0; i < m; i++)
	{
		for (IT j=csrRowPtr[i]; j < csrRowPtr[i+1]; j++)
		{
			dloc[j] = std::atomic_fetch_add(&(atomicColPtr[csrColIdx[j]]), 1);
		}
	}
	cscColPtr[0] = 0;
	for (IT i=0; i <= n; i++)
	{
		cscColPtr[i+1] = atomicColPtr[i] + cscColPtr[i];
	}
	delete [] atomicColPtr;
	
	#pragma omp parallel for schedule(dynamic)
	for (IT i=0; i < m; i++)
	{
		for (IT j=csrRowPtr[i]; j < csrRowPtr[i+1]; j++)
		{
			IT loc = cscColPtr[csrColIdx[j]] + dloc[j];	// dloc[j] (as opposed to a O(N) work array) is what makes this loop race free
			cscRowIdx[loc] = i;
			cscVal[loc] = csrVal[j];
		}
	}
	delete[] dloc;
}
#endif
