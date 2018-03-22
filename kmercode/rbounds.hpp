#ifndef RBOUNDS_HPP
#define RBOUNDS_HPP

#define MIN_PROB 0.001 // 0.1% of unique k-mer

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
int rbounds(int depth, double erate, int klen);

#endif
