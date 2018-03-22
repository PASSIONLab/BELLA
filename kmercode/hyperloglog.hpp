#ifndef HYPERLOGLOG_H
#define HYPERLOGLOG_H

/**
 * @file hyperloglog.hpp
 * @brief HyperLogLog cardinality estimator
 * @date Created 2013/3/20
 * @author Hideaki Ohno
 * Modified by Aydin Buluc to use 64-bit hashes
 */

#include<vector>
#include<cmath>
#include "hash_funcs.h"

#define HASHBITS 64
#define SIGNBIT	0x8000000000000000		// for 32-bit:

static const double pow_2_32 =  4294967296.0;
static const double neg_pow_2_32 = -4294967296.0;

/** @class HyperLogLog
 *  @brief Implement of 'HyperLogLog' estimate cardinality algorithm
 */
class HyperLogLog {
private:
    uint8_t  b;     /// register bit width
    uint32_t m;     /// register size
    double alphaMM; /// alpha * m^2

    uint8_t rho(uint64_t x, uint8_t locb) {
        uint8_t v = 1;
        while (v <= locb && !(x & SIGNBIT)) {
            v++;
            x <<= 1;
        }
        return v;
    }

public:

    std::vector<uint8_t> M;

    /**
     * Constructor
     *
     * @param[in] b_ bit width (register size will be pow(2,b_))
     */
    HyperLogLog(uint8_t b_) : b(b_), m(1<<b), M(m+1,0) {
       double alpha;
       switch(m){
           case 16:
               alpha = 0.673;
               break;
           case 32:
               alpha = 0.697;
               break;
           case 64:
               alpha = 0.709;
               break;
           default:
               alpha = 0.7213/(1.0 + 1.079/m);
       }
       alphaMM = alpha * m * m;
    }

    /**
     * Add element to the estimator
     *
     * @param[in] str string to add
     * @param[in] len length of string
     */
    void add(const char* str, uint32_t len){
        uint64_t hash = MurmurHash3_x64_64(str,len);
        uint32_t index = hash >> ( HASHBITS - b );	// ABAB: use the first b bits to locate the bucket
        uint8_t rank = rho((hash << b), HASHBITS - b);	// ABAB: use the last HASHBITS-b bits to count leading zeros
        if( rank > M[index] ){
            M[index] = rank;
        }
    }

    /**
     * Estimate cardinality value.
     *
     * @return Estimated cardinality value.
     */
    double estimate(){
        double est;
        uint8_t rank = 0;
        double sum = 0.0;
        for (uint32_t i = 0; i < m; i++) {
            sum += 1.0/pow(2.0, M[i]);
        }
        est = alphaMM/sum; // E in the original paper
        if( est <= 2.5 * m ) {
            uint32_t zeros = 0;
            for (uint32_t i = 0; i < m; i++) {
                if (M[i] == 0) {
                    zeros++;
                }
            }
            if( zeros != 0 ) {
                est = m * log((double)m/zeros);
            }
        }
	// ABAB: with 64-bit hash values, large value correction is not needed 
	return est;
    }

};

#endif

