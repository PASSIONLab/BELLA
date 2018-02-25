#ifndef RBOUNDS_HPP
#define RBOUNDS_HPP

#define MIN_PROB 0.001 // 0.1% of unique k-mer

//
// Binomial coefficient function 
//
int bincoef(int n, int k);

//
// Function to select the reliable upper bound
//
int rbounds(int depth, double erate, int klen);

#endif