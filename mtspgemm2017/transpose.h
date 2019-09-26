Function csr2csc_atomic(m, n, nnz, csrRowPtr, csrColIdx, csrVal,
cscColPtr, cscRowIdx, cscVal)
// construct an array of size nnz to record the relative
position of a nonzero element in corresponding column
1 ∗dloc = new int[nnz]();
2 #pragma omp parallel for schedule(dynamic)
3 for i ←0; i < m; i++ do
4 for j ←csrRowPtr[i]; j <csrRowPtr[i+1]; j++ do
5 dloc[j] =
f etch and add(&(cscColPtr[csrColIdx[j] + 1]), 1);
6 pref ix sum(cscColP tr, n + 1);
7 #pragma omp parallel for schedule(dynamic)
8 for i ←0; i < m; i++ do
9 for j ←csrRowPtr[i]; j <csrRowPtr[i+1]; j++ do
10 loc = cscColPtr[csrColIdx[j]] + dloc[j];
11 cscRowIdx[loc] = i;
12 cscVal[loc] = csrV al[j];
13 // sort cscRowIdx,cscVal in each column
14 #pragma omp parallel for schedule(dynamic)
15 for i ←0; i < n; i++ do
16 begin=cscColPtr[i];
17 end =cscColPtr[i + 1];
18 sort key value(being, end, cscRowIdx, cscVal);
19 delete[] dloc;
20 return;
