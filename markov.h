#ifndef BELLA_MARKOV_H_
#define BELLA_MARKOV_H_

#include<bits/stdc++.h>
#include<cmath>
#include<ctgmath>
#include<vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <array>

typedef struct std::vector<std::vector<float>> matrix;

//
//	Expected value of out Markov chain
//

//	adding the entries in the top row, we get the expected number of steps
int getresult(matrix& inverse, const int& dim)
{
	const int row = 0;
	float sum = 0.0;

	for(int col = 0; col < dim; col++)
		sum = sum + inverse[row][col];

	return std::round(sum);
}

void cofactor(matrix& N, matrix& tmp, const int& p, const int& q, const int& dim)
{
	int i = 0, j = 0;
	for(int row = 0; row < dim; row++)
		for(int col = 0; col < dim; col++)
			if(row != p)
				if(col != q)
				{
					tmp[i][j++] = N[row][col];
					// row is filled, increase row index and reset col index
					if(j == dim -1)
					{
						j = 0; i++;
					}
				}
}

//	recursive function for finding determinant of matrix
int determinant(matrix& N, const int& dim) 
{
	float det = 0;

	//	base case: matrix contains single element
	if(dim == 1) 
		return N[0][0]; 

	//	store cofactors
	matrix tmp;
	tmp = std::vector<std::vector<float>>(dim, std::vector<float>(dim, 0.0));
	//	store sifgn multiplier
	float sign = 1;

	//	iterate for each element of first row
	for (int i = 0; i < dim; i++) 
	{
		//	getting cofactor of N[0][f]
		cofactor(N, tmp, 0, i, dim); 
		det += sign * N[0][i] * determinant(tmp, dim - 1);
		//	terms are to be added with alternate sign
		sign = -sign;
	}

	return determinant;
}

//	function to get adjoint of N[dim][dim] in adj[dim][dim]
void adjoint(matrix& N, matrix& adj, const int& dim)
{
	if(dim == 1)
	{
		adj[0][0] = 1.0;
		return;
	}

	int sign = 1;
	//	store cofactors
	matrix tmp;
	tmp = std::vector<std::vector<float>>(dim, std::vector<float>(dim, 0.0));

	for(int i = 0; i < dim; i++)
		for(int j = 0; j < dim; j++)
		{
			cofactor(N, tmp, i, j, dim);
			//	sign of adj[j][i] positive if sum of row and column indexes is even
			sign = ((i + j) % 2 == 0) ? 1: -1;
			//	interchanging rows and columns to get the transpose of the cofactor matrix
			adj[j][i] = (sign) * (determinant(tmp, dim-1));
		}
}

//	function to calculate and store inverse, returns false if matrix is singular
bool inverse(matrix& N, matrix& inverse, const int& dim)
{
	//	find determinant of N
	float det = determinant(N, dim);
	if(det == 0) {
		cout << "Singular matrix, can't find its inverse"; 
		return false;
	}

	//	find adjoint
	matrix adj;
	adj = std::vector<std::vector<float>>(dim, std::vector<float>(dim, 0.0));

	adjoint(N, adj, dim);

	//	find inverse using formula "inverse(N) = adj(N)/det(N)"
	for(int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			inverse[i][j] = adj[i][j]/det;

	return true;
}

//	GG: this assume square matrix
matrix substract(matrix& a, matrix& b, const int& dim)
{
	matrix res;

	for(int i = 0; i < dim; i++)
		for(int j = 0; j < dim; j++) {
			res[i][j] = a[i][j] - b[i][j];
	}
	return res;
}

void generaten(matrix& Q, matrix& I, matrix& N, const int& dim) {
	N = substract(I, Q, dim);
}

void generatei(matrix& I, const int& dim)
{
	//	GG: Let Q be the sub-matrix of P without the rows and columns of any absorbing states
	I = std::vector<std::vector<float>>(dim, std::vector<float>(dim, 0.0));
	for(int i = 0; i < dim; i++)
		I[i][i] = 1.0;
}

void generateq(const matrix& P, matrix& Q, const int& dim)
{
	//	GG: Let Q be the sub-matrix of P without the rows and columns of any absorbing states
	Q = std::vector<std::vector<float>>(dim, std::vector<float>(dim, 0.0));
	Q = P;
	Q.resize(dim - 1, std::vector<float>(dim - 1));
}

int markovstep(const float& probability, const int& kmersize)
{
	const int dim  = kmersize + 1;
	int res = -1;

	//	GG: transition matrix
	matrix P, Q, I, N, inversemat;
	P = std::vector<std::vector<float>>(dim, std::vector<float>(dim, 0.0));

	for(int i = 0; i < kmersize; i++)
	{
		P[i][0]   = (1 - std::pow(probability, 2));
		P[i][i+1] = std::pow(probability, 2);
	}
	P[kmersize][kmersize] = 1.0;

	//	GG: Let Q be the sub-matrix of P without the rows and columns of any absorbing states
	generateq(P, Q, dim);
	//	GG: generate identity matrix
	generatei(I, kmersize);
	//	GG: generate fundamental matrix
	generaten(Q, I, N, kmersize);

	if(inverse(N, inversemat, kmersize))
		res = getresult(inversemat, kmersize);

	return res;	//	expected overlap length to get a correct kmer
}

#endif