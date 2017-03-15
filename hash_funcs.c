#include <stdint.h>

#include "hash_funcs.h"

uint32_t rotl32(uint32_t x, int8_t r)
{
	return (x << r) | (x >> (32 - r));
}

uint64_t rotl64(uint64_t x, int8_t r)
{
	return (x << r) | (x >> (64 - r));
}

#define ROTL32(x,y) rotl32(x,y)
#define ROTL64(x,y) rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

#define getblock(p, i) (p[i])

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

uint64_t fmix64(uint64_t k)
{
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xff51afd7ed558ccd);
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
	k ^= k >> 33;

	return k;
}

//-----------------------------------------------------------------------------

void MurmurHash3_x64_128(const void * key, const uint32_t len, const uint32_t seed, void * out)
{
	const uint8_t * data = (const uint8_t*)key;
	const uint32_t nblocks = len / 16;
	int32_t i;

	uint64_t h1 = seed;
	uint64_t h2 = seed;

	uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
	uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

	//----------
	// body

	const uint64_t * blocks = (const uint64_t *)(data);

	for(i = 0; i < nblocks; i++)
	{
		uint64_t k1 = getblock(blocks,i*2+0);
		uint64_t k2 = getblock(blocks,i*2+1);

		k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;

		h1 = ROTL64(h1,27); h1 += h2; h1 = h1*5+0x52dce729;

		k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

		h2 = ROTL64(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
	}

	//----------
	// tail

	const uint8_t * tail = (const uint8_t*)(data + nblocks*16);

	uint64_t k1 = 0;
	uint64_t k2 = 0;

	switch(len & 15)
	{
	case 15: k2 ^= (uint64_t)(tail[14]) << 48;
	case 14: k2 ^= (uint64_t)(tail[13]) << 40;
	case 13: k2 ^= (uint64_t)(tail[12]) << 32;
	case 12: k2 ^= (uint64_t)(tail[11]) << 24;
	case 11: k2 ^= (uint64_t)(tail[10]) << 16;
	case 10: k2 ^= (uint64_t)(tail[ 9]) << 8;
	case  9: k2 ^= (uint64_t)(tail[ 8]) << 0;
		k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

	case  8: k1 ^= (uint64_t)(tail[ 7]) << 56;
	case  7: k1 ^= (uint64_t)(tail[ 6]) << 48;
	case  6: k1 ^= (uint64_t)(tail[ 5]) << 40;
	case  5: k1 ^= (uint64_t)(tail[ 4]) << 32;
	case  4: k1 ^= (uint64_t)(tail[ 3]) << 24;
	case  3: k1 ^= (uint64_t)(tail[ 2]) << 16;
	case  2: k1 ^= (uint64_t)(tail[ 1]) << 8;
	case  1: k1 ^= (uint64_t)(tail[ 0]) << 0;
		k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
	};

	//----------
	// finalization

	h1 ^= len; h2 ^= len;

	h1 += h2;
	h2 += h1;

	h1 = fmix64(h1);
	h2 = fmix64(h2);

	h1 += h2;
	h2 += h1;

	((uint64_t*)out)[0] = h1;
	((uint64_t*)out)[1] = h2;
}


//-----------------------------------------------------------------------------
// If we need a smaller hash value, it's faster to just use a portion of the 
// 128-bit hash

uint32_t MurmurHash3_x64_32(const void * key, uint32_t len)
{
	uint32_t temp[4];
	MurmurHash3_x64_128(key, len, 313, temp);
	return temp[0];
}

//----------

uint64_t MurmurHash3_x64_64(const void * key, uint32_t len)
{
	uint64_t temp[2];
	MurmurHash3_x64_128(key, len, 313, temp);
	return temp[0];
} 

//-----------------------------------------------------------------------------

// http://murmurhash.googlepages.com/MurmurHash2.cpp
uint64_t murmur_hash2_32(const char * key, uint32_t len)
{
	// 'm' and 'r' are mixing constants generated offline.
	// They're not really 'magic', they just happen to work well.
	const uint32_t m = 0x5bd1e995;
	const int32_t r = 24;
	// Initialize the hash to a 'random' value
	uint32_t seed = 0x3FB0BB5F;
	uint32_t h = seed ^ len;
	// Mix 4 bytes at a time into the hash
	const uint8_t * data = (const uint8_t *)key;
	while(len >= 4)	{
		uint32_t k = *(uint32_t *)data;
		k *= m; 
		k ^= k >> r; 
		k *= m; 
		
		h *= m; 
		h ^= k;

		data += 4;
		len -= 4;
	}
	// Handle the last few bytes of the input array
	switch(len)	{
	case 3: h ^= data[2] << 16;
	case 2: h ^= data[1] << 8;
	case 1: h ^= data[0];
	        h *= m;
	};
	// Do a few final mixes of the hash to ensure the last few
	// bytes are well-incorporated.
	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;
	return h;
} 

uint64_t murmur_hash2_64(const void * key, uint32_t len)
{
	uint32_t seed = 0x3FB0BB5F;
	const uint64_t m = BIG_CONSTANT(0xc6a4a7935bd1e995);
	const int32_t r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while(data != end)
	{
		uint64_t k = *data++;

		k *= m; 
		k ^= k >> r; 
		k *= m; 
    
		h ^= k;
		h *= m; 
	}

	const uint8_t * data2 = (const uint8_t*)data;

	switch(len & 7)
	{
	case 7: h ^= ((uint64_t) data2[6]) << 48;
	case 6: h ^= ((uint64_t) data2[5]) << 40;
	case 5: h ^= ((uint64_t) data2[4]) << 32;
	case 4: h ^= ((uint64_t) data2[3]) << 24;
	case 3: h ^= ((uint64_t) data2[2]) << 16;
	case 2: h ^= ((uint64_t) data2[1]) << 8;
	case 1: h ^= ((uint64_t) data2[0]);
		h *= m;
	};
 
	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
} 

#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)\
                       +(uint32_t)(((const uint8_t *)(d))[0]) )
#endif

uint32_t SuperFastHash (const char * data, int len) {
    uint32_t hash = len, tmp;
    int rem;

    if (len <= 0 || !data) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (data);
                hash ^= hash << 16;
                hash ^= data[sizeof (uint16_t)] << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (data);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += *data;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}


