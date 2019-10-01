#include <atomic>

/** 
 ** while the name suggests CSR -> CSC, 
 ** one can simply use this to transpose a matrix in CSR format to its transpose also in CSR format
 ** or vice versa: transpose a matrix in CSC format to its transpose also in CSC format
 ** just read CSR related params as the "old matrix" and CSC related params as the "output/new matrix"
 **/
template <class IT, class NT>
void csr2csc_atomic (IT m, IT n, IT nnz, IT * csrRowPtr, IT * csrColIdx, NT * csrVal,
			IT * cscColPtr, IT * cscRowIdx, NT * cscVal)
{
	std::atomic<int> * atomicColPtr = new std::atomic<int>[n+1];
	for (IT i=0; i < m; i++)
	{
		atomic_init(atomicColPtr[i], 0);
	}

	// construct an array of size nnz to record the relative
	// position of a nonzero element in corresponding column
	IT * dloc = new IT[nnz]();
	#pragma omp parallel for schedule(dynamic)
	for (IT i=0; i < m; i++)
	{
		for (IT j=csrRowPtr[i]; j < csrRowPtr[i+1]; j++)
		{
			dloc[j] = std::atomic_fetch_add(&(atomicColPtr[csrColIdx[j]+1]), 1);
		}
	}
	// The output iterator result is allowed to be the same iterator as the input iterator first,
	// so that partial sums may be computed in place; but note we are writing to the non-atomic
	std::partial_sum (atomicColPtr, atomicColPtr+n+1, cscColPtr);
	delete [] atomicColPtr;
	
	#pragma omp parallel for schedule(dynamic)
	for (IT i=0; i < m; i++)
	{
		for (IT j=csrRowPtr[i]; j < csrRowPtr[i+1]; j++)
		{
			IT loc = cscColPtr[csrColIdx[j]] + dloc[j];
			cscRowIdx[loc] = i;
			cscVal[loc] = csrVal[j];
		}
	}
	delete[] dloc;

	// sort cscRowIdx,cscVal in each column
	#pragma omp parallel for schedule(dynamic)
	for (IT i=0; i < n; i++)
	{
		// can we do this without using std::pair?
	}
}
