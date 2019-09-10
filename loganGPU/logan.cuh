//==================================================================
// Title:  Cuda x-drop seed-and-extend alignment algorithm
// Author: G. Guidi, A. Zeni
// Date:   6 March 2019
//==================================================================

#ifndef __LOGAN_CUH__
#define __LOGAN_CUH__

#include<vector>
#include<iostream>
#include<chrono>
#include<numeric>
#include<functional>
#include<iterator>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include"seed.cuh"
#include"score.cuh"
#include"functions.cuh"
//#include <boost/algorithm/string.hpp>
//#include <cub/block/block_load.cuh>
//#include <cub/block/block_store.cuh>
//#include <cub/block/block_reduce.cuh>
//#include <cub/cub.cuh>
//#include<boost/array.hpp>

// __inline__ __device__ short simple_max(int *antidiag, int &dim, int &offset);
// __inline__ __device__ void warpReduce(volatile int *input, int myTId);
// __inline__ __device__ short warpReduceMax(short val);
// __inline__ __device__ short blockReduceMax(short val);
// __inline__ __device__ int reduce_max(int *antidiag);
// __inline__ __device__ void updateExtendedSeedL(SeedL& seed, ExtensionDirectionL direction, int &cols, int &rows, int &lowerDiag, int &upperDiag);
// __inline__ __device__ void computeAntidiag(int *antiDiag1, int *antiDiag2, int *antiDiag3, char* querySeg, char* databaseSeg, int &best, int &scoreDropOff, int &cols, int &rows, int &minCol, int &maxCol, int &antiDiagNo, int &offset1, int &offset2, ExtensionDirectionL direction);
// __inline__ __device__ void calcExtendedLowerDiag(int *lowerDiag, int const &minCol, int const &antiDiagNo);
// __inline__ __device__ void calcExtendedUpperDiag(int *upperDiag, int const &maxCol, int const &antiDiagNo);
// __inline__ __device__ void initAntiDiag3(int *antiDiag3, int *a3size, int const &offset, int const &maxCol, int const &antiDiagNo, int const &minScore, int const &gapCost, int const &undefined);
// __inline__ __device__ void initAntiDiags(int *antiDiag1,int *antiDiag2,int *antiDiag3,int *a2size,int *a3size,int const &dropOff,int const &gapCost,int const &undefined);
// __global__ void extendSeedLGappedXDropOneDirection( SeedL *seed, char *querySegArray, char *databaseSegArray, ExtensionDirectionL direction, int scoreDropOff, int *res, int *qL, int *dbL, int *offsetQuery, int *offsetTarget);
// inline void extendSeedL(vector<SeedL> &seeds, ExtensionDirectionL direction, vector<string> &target, vector<string> &query, vector<ScoringSchemeL> &penalties, int const& XDrop, int const& kmer_length);

#endif
