#ifndef BELLA_KMERCODE_BOUND_H_
#define BELLA_KMERCODE_BOUND_H_

//	GG:	select upper and lower reliable bounds
int computeUpper(int myCoverage, double errorRate, int kmerSize, double minProbability);
int computeLower(int myCoverage, double errorRate, int kmerSize, double minProbability);

#endif
