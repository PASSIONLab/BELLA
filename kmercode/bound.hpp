#ifndef BELLA_KMERCODE_BOUND_H_
#define BELLA_KMERCODE_BOUND_H_

#define MIN_PROB 0.002 // 0.2% cumulative sum

/**
 * @brief bincoef computes the binomial coefficient
 * @param n
 * @param k
 * @return
 */
int bincoef(int n, int k);

/**
 * @brief rbounds selects the reliable upper bound
 * @param depth
 * @param erate
 * @param klen
 * @return
 */
int computeUpper(int depth, double erate, int klen);
int computeLower(int depth, double erate, int klen);

#endif
