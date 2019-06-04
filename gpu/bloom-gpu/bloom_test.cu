#include <stdio.h>

#include <cstdint>
#include <cstdlib>
#include <chrono>
#include <iostream>

typedef std::chrono::high_resolution_clock Clock;

#include "nvbio/bloom_filter.h"
#include "nvbio/types.h"

#include "fhash.h"

// populate a filter in parallel
template <typename bloom_filter_type, typename T>
__global__ 
void 
populate_kernel(nvbio::uint32 N, T* arr, bloom_filter_type filter)
{
  nvbio::uint32 i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < N)
    filter.insert(arr[i]);
}


int
main (int argc, char **argv)
{
	srand(7);
	int32_t N = atoi(argv[1]);

	uint32_t *h_data = (uint32_t *) malloc(sizeof(*h_data) * N);
	for (nvbio::uint32 i = 0; i < N; ++i)
	{
		h_data[i] = rand();
		// std::cout << "elem " << i << " " << h_data[i].n << "\n";
	}


	auto t1 = Clock::now();
	uint32_t *d_data = NULL;
	cudaMalloc((void **)&d_data, N * sizeof(*d_data));
	cudaMemcpy(d_data, h_data, N * sizeof(*d_data), cudaMemcpyHostToDevice);
	
	typedef nvbio::bloom_filter<4, RSHash<uint32_t>, ElfHash<uint32_t>, nvbio::uint32*>
		bloom_filter_type;

	int nfilter_elems = N / 4;
	uint32_t *d_filter_storage = NULL;
	cudaMalloc((void **)&d_filter_storage,
			   nfilter_elems * sizeof(*d_filter_storage));
	cudaMemset(d_filter_storage, 0,
			   nfilter_elems * sizeof(*d_filter_storage));
	
	bloom_filter_type d_filter(nfilter_elems * 32, d_filter_storage);
	
	int nblocks = (N + 1024) / 1024;
	populate_kernel<<<nblocks,1024>>>(N, d_data, d_filter);
	cudaDeviceSynchronize();

	auto t2 = Clock::now();
	double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
	duration = duration / 1e6;
	
	printf("bloom filter construction took %.2f milliseconds for %d elems\n",
		   duration, N);
	
	uint32_t *h_filter_storage = (uint32_t *)
		malloc(sizeof(*h_filter_storage) * nfilter_elems);
	cudaMemcpy(h_filter_storage, d_filter_storage,
			   nfilter_elems * sizeof(*h_filter_storage),
			   cudaMemcpyDeviceToHost);
	bloom_filter_type h_filter(nfilter_elems * 32, h_filter_storage);



	// std::cout << "Bloom filter on CPU\n" << std::flush;
	// memset(h_filter_storage, nfilter_elems, sizeof(*h_filter_storage));
	// bloom_filter_type h_filter2(nfilter_elems * 32, h_filter_storage);
	// for (int i = 0; i < N; ++i)
	// 	h_filter2.insert(h_data[i]);


	////////////////////////////////////////////////////////////////////////////
	// uint32_t kex;
	// kex.n = 1045618677;
	// kex.len = 4;
	// std::cout << h_filter.has(kex) << "\n";
	
	// kex.n = 8989;
	// kex.len = 4;
	// std::cout << h_filter.has(kex) << "\n";

	// kex.n = 781032961;
	// kex.len = 4;
	// std::cout << h_filter.has(kex) << "\n";
	////////////////////////////////////////////////////////////////////////////


}
