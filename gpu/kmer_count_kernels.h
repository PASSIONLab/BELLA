#include <cstdint>

#include "kmercode/Kmer.hpp"

template<typename T>
struct minus_ptr
{
	__host__ __device__
	ptrdiff_t
	operator()(T a, T b)
	{
		return a - b;
	}
};

struct equal_kmer
{
  	__host__ __device__
	bool
	operator()(uint64_t *a, uint64_t *b)
	{
		bool r = true;
		for (size_t i = 0; i < N_LONGS; ++i)
		{
			r &= *a == *b;
			++a;
			++b;
		}
		return r;
	}
};

struct compare_kmer
{
  	__host__ __device__
	bool
	operator()(uint64_t *a, uint64_t *b)
	{
		bool r = false;
		for (size_t i = 0; i < N_LONGS; ++i)
		{
			if (*a < *b)
				return true;
			if (*a > *b)
				return false;
			++a;
			++b;
		}
		return false;
	}
};


struct printf_functor
{
	 __host__ __device__
  	void operator()(int x)
	{
		printf("%d\n", x);
	}
};


// populate a filter in parallel
template <typename bloom_filter_type, typename T>
__global__ 
void 
populate_kernel(uint32_t N, T* arr, bloom_filter_type filter,
				uint8_t *d_kmer_pass, uint64_t **d_kmer_ptrs)
{
  	uint32_t i = threadIdx.x + blockIdx.x * blockDim.x;	
	if (i < N)
	{
		uint32_t idx = i * N_LONGS;
		bool isinfilter = filter.has(&(arr[idx]));
		// printf("%d %d is in filter %d\n", i, N, isinfilter);
		if (isinfilter)
		{
			// printf("thread %d already have it...\n", i);
			d_kmer_pass[i] = 1;
			d_kmer_ptrs[i]  = &(arr[idx]);
		}
		else
		{
			// printf("thread %d inserting...\n", i);
			filter.insert(&(arr[idx]));
		}
	}
}
