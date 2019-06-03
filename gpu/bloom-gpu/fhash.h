// WARNING: Do not use for anything else but kmers

#include "kmercode/Kmer.hpp"


template<typename T>
struct ElfHash
{
	__host__ __device__
    unsigned int
	operator()(T k)
    {
		// printf("in ElfHash\n");
		char *ptr = (char *)k;

		unsigned int hash = 0;
		unsigned int x    = 0;
		unsigned int i    = 0;

		for (i = 0; i < N_BYTES; ++ptr, ++i)
		{
			// printf("%d ", *ptr);
			hash = (hash << 4) + (*ptr);
			if ((x = hash & 0xF0000000L) != 0)
			{
				hash ^= (x >> 24);
			}
			hash &= ~x;
		}
		// printf("hash: %u\n", hash);
		return hash;
	}

	__host__ __device__
    unsigned int
	operator()(T const k) const
    {
		// printf("in ElfHash\n");
		const char *ptr = (char *)k;

		unsigned int hash = 0;
		unsigned int x    = 0;
		unsigned int i    = 0;

		for (i = 0; i < N_BYTES; ++ptr, ++i)
		{
			hash = (hash << 4) + (*ptr);
			if ((x = hash & 0xF0000000L) != 0)
			{
				hash ^= (x >> 24);
			}
			hash &= ~x;
		}
		// printf("hash: %u\n", hash);
		return hash;
	}
};

template<typename T>
struct RSHash
{
	// @TODO should ref
    __host__ __device__
    unsigned int
	operator()(T k)
    {
		// printf("in RSHash\n");
		char *ptr = (char *)k;
		
		unsigned int	 b    = 378551;
		unsigned int	 a    = 63689;
		unsigned int	 hash = 0;
		unsigned int	 i    = 0;		
		for (i = 0; i < N_BYTES; ++ptr, ++i)
		{
			hash = hash * a + (*ptr);
			a    = a * b;
		}
		// printf("hash: %u\n", hash);
		return hash;
	}

	__host__ __device__
    unsigned int
	operator()(T const k) const
    {
		// printf("in RSHash\n");
		const char *ptr = (char *)k;
		
		unsigned int	 b    = 378551;
		unsigned int	 a    = 63689;
		unsigned int	 hash = 0;
		unsigned int	 i    = 0;		
		for (i = 0; i < N_BYTES; ++ptr, ++i)
		{
			hash = hash * a + (*ptr);
			a    = a * b;
		}
		// printf("hash: %u\n", hash);
		return hash;
	}	
};
