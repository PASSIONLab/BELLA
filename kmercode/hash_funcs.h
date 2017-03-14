#ifndef HASH_FUNCS_H
#define HASH_FUNCS_H

#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif


uint32_t rotl32(uint32_t x, int8_t r);

uint64_t rotl64(uint64_t x, int8_t r);

#define ROTL32(x,y) rotl32(x,y)
#define ROTL64(x,y) rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

#define getblock(p, i) (p[i])

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

uint64_t fmix64(uint64_t k);

//-----------------------------------------------------------------------------

void MurmurHash3_x64_128(const void * key, const uint32_t len, const uint32_t seed, void * out);


//-----------------------------------------------------------------------------
// If we need a smaller hash value, it's faster to just use a portion of the 
// 128-bit hash

uint32_t MurmurHash3_x64_32(const void * key, uint32_t len);

//----------

uint64_t MurmurHash3_x64_64(const void * key, uint32_t len);

//-----------------------------------------------------------------------------

// http://murmurhash.googlepages.com/MurmurHash2.cpp
uint64_t murmur_hash2_32(const char * key, uint32_t len);

uint64_t murmur_hash2_64(const void * key, uint32_t len);

uint32_t SuperFastHash (const char * data, int len);

#ifdef __cplusplus
}
#endif



#endif // HASH_FUNCS_H
