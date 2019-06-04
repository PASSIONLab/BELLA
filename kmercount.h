#ifndef BELLA_KMERCOUNT_H_
#define BELLA_KMERCOUNT_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <array>
#include <tuple>
#include <queue>
#include <memory>
#include <stack>
#include <numeric>
#include <functional>
#include <cstring>
#include <string.h>
#include <math.h>
#include <cassert>
#include <ios>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <map>
#include <unordered_map>
#include <omp.h>
#include <chrono>

#include "libcuckoo/cuckoohash_map.hh"
#include "libbloom/bloom64.h"

#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include "kmercode/Buffer.h"
#include "kmercode/common.h"
#include "kmercode/fq_reader.h"
#include "kmercode/ParallelFASTQ.h"
#include "kmercode/bound.hpp"
#include "kmercode/hyperloglog.hpp"
#include "mtspgemm2017/common.h"

#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <cuda_profiler_api.h>
#include "gpu/kmer_count_kernels.h"
#include "gpu/bloom-gpu/fhash.h"
#include "gpu/bloom-gpu/nvbio/bloom_filter.h"
#include "gpu/bloom-gpu/nvbio/types.h"
typedef std::chrono::high_resolution_clock Clock;


using namespace std;
#define ASCIIBASE 33 // Pacbio quality score ASCII BASE
#ifndef PRINT
#define PRINT
#endif

typedef cuckoohash_map<Kmer, int> dictionary_t; // <k-mer && reverse-complement, #kmers>

struct filedata {

    char filename[MAX_FILE_PATH];
    size_t filesize;
};

/**
 * @brief GetFiles
 * @param filename
 * @return
 */
vector<filedata>  GetFiles(char *filename) {
    int64_t totalsize = 0;
    int numfiles = 0;
    std::vector<filedata> filesview;
    
    filedata fdata;
    ifstream allfiles(filename);
    if(!allfiles.is_open()) {
        cerr << "Could not open " << filename << endl;
        exit(1);
    }
    allfiles.getline(fdata.filename,MAX_FILE_PATH);
    while(!allfiles.eof())
    {
        struct stat st;
        stat(fdata.filename, &st);
        fdata.filesize = st.st_size;
        
        filesview.push_back(fdata);
        cout << filesview.back().filename << ": " << filesview.back().filesize / (1024*1024) << " MB" << endl;
        allfiles.getline(fdata.filename,MAX_FILE_PATH);
        totalsize += fdata.filesize;
        numfiles++;
    }
    return filesview;
}

/**
 * @brief JellyFishCount
 * @param kmer_file
 * @param countsreliable_jelly
 * @param lower
 * @param upper
 */
void JellyFishCount(char *kmer_file, dictionary_t & countsreliable_jelly, int lower, int upper) 
{
    ifstream filein(kmer_file);
    string line;
    int elem;
    string kmerstr;    
    Kmer kmerfromstr;
    
    // double kdict = omp_get_wtime();
    // Jellyfish file contains all the k-mers from fastq(s)
    // It is not filtered beforehand
    // A k-mer and its reverse complement are counted separately
    dictionary_t countsjelly;
    if(filein.is_open()) 
    { 
            while(getline(filein, line)) {
                if(line.length() == 0)
                    break;

                string substring = line.substr(1);
                elem = stoi(substring);
                getline(filein, kmerstr);   
                //kmerfromstr.set_kmer(kmerstr.c_str());

                auto updatecountjelly = [&elem](int &num) { num+=elem; };
                // If the number is already in the table, it will increment its count by the occurrence of the new element. 
                // Otherwise it will insert a new entry in the table with the corresponding k-mer occurrence.
                countsjelly.upsert(kmerfromstr.rep(), updatecountjelly, elem);      
            }
    } else std::cout << "Unable to open the input file\n";
    filein.close();
    //cout << "jellyfish file parsing took: " << omp_get_wtime()-kdict << "s" << endl;

    // Reliable k-mer filter on countsjelly
    int kmer_id = 0;
    auto lt = countsjelly.lock_table(); // our counting
    for (const auto &it : lt) 
        if (it.second >= lower && it.second <= upper)
        {
            countsreliable_jelly.insert(it.first,kmer_id);
            ++kmer_id;
        }
    lt.unlock(); // unlock the table
    // Print some information about the table
    cout << "Entries within reliable range Jellyfish: " << countsreliable_jelly.size() << std::endl;    
    //cout << "Bucket count Jellyfish: " << countsjelly.bucket_count() << std::endl;
    //cout << "Load factor Jellyfish: " << countsjelly.load_factor() << std::endl;
    countsjelly.clear(); // free 
}

/**
 * @brief DeNovoCount
 * @param allfiles
 * @param countsreliable_denovo
 * @param lower
 * @param upper
 * @param kmer_len
 * @param upperlimit
 */
void DeNovoCount_cpu(vector<filedata> & allfiles, dictionary_t & countsreliable_denovo, int & lower, int & upper, int kmer_len, int depth, double & erate, size_t upperlimit /* memory limit */, BELLApars & b_parameters)
{
    vector < vector<Kmer> > allkmers(MAXTHREADS);
    vector < vector<double> > allquals(MAXTHREADS);
    vector < HyperLogLog > hlls(MAXTHREADS, HyperLogLog(12));   // std::vector fill constructor

    double denovocount = omp_get_wtime();
    double cardinality;
    size_t totreads = 0;

    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) 
    {
        #pragma omp parallel
        {
            ParallelFASTQ *pfq = new ParallelFASTQ();
            pfq->open(itr->filename, false, itr->filesize);

            vector<string> seqs;
            vector<string> quals;
            vector<string> nametags;
            size_t tlreads = 0; // thread local reads

            size_t fillstatus = 1;
            while(fillstatus) 
            { 
                fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
                size_t nreads = seqs.size();

                for(int i=0; i<nreads; i++) 
                {
                    // remember that the last valid position is length()-1
                    int len = seqs[i].length();
                    double rerror = 0.0;

                    for(int j=0; j<=len-kmer_len; j++)  
                    {
                        std::string kmerstrfromfastq = seqs[i].substr(j, kmer_len);
                        Kmer mykmer(kmerstrfromfastq.c_str());
                        Kmer lexsmall = mykmer.rep();
                        allkmers[MYTHREAD].push_back(lexsmall);
                        hlls[MYTHREAD].add((const char*) lexsmall.getBytes(), lexsmall.getNumBytes());

            		if(b_parameters.skipEstimate == false)
            		{
                        	// accuracy
                       		int bqual = (int)quals[i][j] - ASCIIBASE;
                        	double berror = pow(10,-(double)bqual/10);
                        	rerror += berror;
            		}

                    }
		    if(b_parameters.skipEstimate == false)
		    {
                    	// remaining k qual position accuracy
                    	for(int j=len-kmer_len+1; j < len; j++)
                    	{
                        	int bqual = (int)quals[i][j] - ASCIIBASE;
                        	double berror = pow(10,-(double)bqual/10);
                        	rerror += berror;
                    	}
                    	rerror = rerror / len;
                    	allquals[MYTHREAD].push_back(rerror);
		    }
                } // for(int i=0; i<nreads; i++)
                tlreads += nreads;
            } //while(fillstatus) 
            delete pfq;

            #pragma omp critical
            totreads += tlreads;
        }
    }

    // Error estimation
    if(b_parameters.skipEstimate == false)
    {
        erate = 0.0; // reset to 0 here, otherwise it cointains default or user-defined values
        #pragma omp for reduction(+:erate)
        for (int i = 0; i < MAXTHREADS; i++) 
            {
                double temp = std::accumulate(allquals[i].begin(),allquals[i].end(), 0.0);
                erate += temp/(double)allquals[i].size();
            }
        erate = erate / (double)MAXTHREADS;
    }

    // HLL reduction (serial for now) to avoid double iteration
    for (int i = 1; i < MAXTHREADS; i++) 
    {
        std::transform(hlls[0].M.begin(), hlls[0].M.end(), hlls[i].M.begin(), hlls[0].M.begin(), [](uint8_t c1, uint8_t c2) -> uint8_t{ return std::max(c1, c2); });
    }
    cardinality = hlls[0].estimate();

    double load2kmers = omp_get_wtime(); 
    cout << "Initial parsing, error estimation, and k-mer loading took: " << load2kmers - denovocount << "s\n" << endl;

    const double desired_probability_of_false_positive = 0.05;
    struct bloom * bm = (struct bloom*) malloc(sizeof(struct bloom));
    bloom_init64(bm, cardinality * 1.1, desired_probability_of_false_positive);

#ifdef PRINT
    cout << "Cardinality estimate is " << cardinality << endl;
    cout << "Table size is: " << bm->bits << " bits, " << ((double)bm->bits)/8/1024/1024 << " MB" << endl;
    cout << "Optimal number of hash functions is: " << bm->hashes << endl;
#endif

    dictionary_t countsdenovo;


	uint64_t totkmers = 0;
	for (int i = 0; i < MAXTHREADS; ++i)
		totkmers += allkmers[i].size();

	auto t1 = Clock::now();
	#pragma omp parallel
    {       
    	for(auto v:allkmers[MYTHREAD])
    	{
        	bool inBloom = (bool) bloom_check_add(bm, v.getBytes(), v.getNumBytes(),1);
        	if(inBloom) countsdenovo.insert(v, 0);
    	}
    }


    double firstpass = omp_get_wtime();
    cout << "First pass of k-mer counting took: " << firstpass - load2kmers << "s" << endl;

    free(bm); // release bloom filter memory

    // in this pass, only use entries that already are in the hash table
    auto updatecount = [](int &num) { ++num; };
#pragma omp parallel
    {       
    	for(auto v:allkmers[MYTHREAD])
    	{
        	// does nothing if the entry doesn't exist in the table
        	countsdenovo.update_fn(v,updatecount);
    	}
    }
    cout << "Second pass of k-mer counting took: " << omp_get_wtime() - firstpass << "s\n" << endl;

		auto t2 = Clock::now();
	double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
	duration = duration / 1e6;
	printf("bloom filter insert/query took %.2f milliseconds on CPU for %d kmers\n",
		   duration, totkmers);

    //cout << "countsdenovo.size() " << countsdenovo.size() << endl;
    // Reliable bounds computation using estimated error rate from phred quality score
    lower = computeLower(depth, erate, kmer_len);
    upper = computeUpper(depth, erate, kmer_len);

    // Reliable k-mer filter on countsdenovo
	uint32_t nkmersintable = 0;
    int kmer_id_denovo = 0;
    auto lt = countsdenovo.lock_table(); // our counting
    for (const auto &it : lt)
	{
		++nkmersintable;
        if (it.second >= lower && it.second <= upper)
        {
            countsreliable_denovo.insert(it.first,kmer_id_denovo);
            ++kmer_id_denovo;
        }
	}
    lt.unlock(); // unlock the table

	cout << "#kmers in table: " << nkmersintable << "\n";

    // Print some information about the table
    if (countsreliable_denovo.size() == 0)
    {
        cout << "BELLA terminated: 0 entries within reliable range (reduce k-mer length)\n" << endl;
        // exit(0);
    } 
    else 
    {
        cout << "Entries within reliable range: " << countsreliable_denovo.size() << endl;
    }
    //cout << "Bucket count: " << countsdenovo.bucket_count() << std::endl;
    //cout << "Load factor: " << countsdenovo.load_factor() << std::endl;
    countsdenovo.clear(); // free

}



/**
 * @brief DeNovoCount
 * @param allfiles
 * @param countsreliable_denovo
 * @param lower
 * @param upper
 * @param kmer_len
 * @param upperlimit
 */
void
DeNovoCount(vector<filedata> & allfiles,
			dictionary_t & countsreliable_denovo,
			int & lower,
			int & upper,
			int kmer_len,
			int depth,
			double & erate,
			size_t upperlimit /* memory limit */,
			BELLApars & b_parameters)
{
    vector < vector<Kmer> > allkmers(MAXTHREADS);
    vector < vector<double> > allquals(MAXTHREADS);
    vector < HyperLogLog > hlls(MAXTHREADS, HyperLogLog(12));   // std::vector fill constructor

    double denovocount = omp_get_wtime();
    double cardinality;
    size_t totreads = 0;

	// cout << "sizeof Kmer class: " << sizeof(Kmer) << "\n" << std::flush;
	
    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) 
    {
        #pragma omp parallel
        {
            ParallelFASTQ *pfq = new ParallelFASTQ();
            pfq->open(itr->filename, false, itr->filesize);

            vector<string> seqs;
            vector<string> quals;
            vector<string> nametags;
            size_t tlreads = 0; // thread local reads

            size_t fillstatus = 1;
            while(fillstatus) 
            { 
                fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
                size_t nreads = seqs.size();

                for(int i=0; i<nreads; i++) 
                {
                    // remember that the last valid position is length()-1
                    int len = seqs[i].length();
                    double rerror = 0.0;

                    for(int j=0; j<=len-kmer_len; j++)  
                    {
                        std::string kmerstrfromfastq = seqs[i].substr(j, kmer_len);
                        Kmer mykmer(kmerstrfromfastq.c_str());
                        Kmer lexsmall = mykmer.rep();
                        allkmers[MYTHREAD].push_back(lexsmall);
                        hlls[MYTHREAD].add((const char*) lexsmall.getBytes(), lexsmall.getNumBytes());

            		if(b_parameters.skipEstimate == false)
            		{
                        	// accuracy
                       		int bqual = (int)quals[i][j] - ASCIIBASE;
                        	double berror = pow(10,-(double)bqual/10);
                        	rerror += berror;
            		}

                    }
		    if(b_parameters.skipEstimate == false)
		    {
                    	// remaining k qual position accuracy
                    	for(int j=len-kmer_len+1; j < len; j++)
                    	{
                        	int bqual = (int)quals[i][j] - ASCIIBASE;
                        	double berror = pow(10,-(double)bqual/10);
                        	rerror += berror;
                    	}
                    	rerror = rerror / len;
                    	allquals[MYTHREAD].push_back(rerror);
		    }
                } // for(int i=0; i<nreads; i++)
                tlreads += nreads;
            } //while(fillstatus) 
            delete pfq;

            #pragma omp critical
            totreads += tlreads;
        }
    }

    // Error estimation
    if(b_parameters.skipEstimate == false)
    {
        erate = 0.0; // reset to 0 here, otherwise it cointains default or user-defined values
        #pragma omp for reduction(+:erate)
        for (int i = 0; i < MAXTHREADS; i++) 
            {
                double temp = std::accumulate(allquals[i].begin(),allquals[i].end(), 0.0);
                erate += temp/(double)allquals[i].size();
            }
        erate = erate / (double)MAXTHREADS;
    }

    // HLL reduction (serial for now) to avoid double iteration
    for (int i = 1; i < MAXTHREADS; i++) 
    {
        std::transform(hlls[0].M.begin(), hlls[0].M.end(), hlls[i].M.begin(), hlls[0].M.begin(), [](uint8_t c1, uint8_t c2) -> uint8_t{ return std::max(c1, c2); });
    }
    cardinality = hlls[0].estimate();

    double load2kmers = omp_get_wtime(); 
    cout << "Initial parsing, error estimation, and k-mer loading took: " << load2kmers - denovocount << "s\n" << endl;


	////////////////////////////////////////////////////////////////////////////
	uint64_t totkmers = 0;
	for (int i = 0; i < MAXTHREADS; ++i)
		totkmers += allkmers[i].size();
	assert (N_LONGS == 4);
	cout << "Copying kmers to gpu...\n";
	uint64_t *h_kmers = (uint64_t *)
		malloc(sizeof(*h_kmers) * totkmers * N_LONGS);
	uint64_t tmp = 0;
	for (int i = 0; i < MAXTHREADS; i++)
	{
		for (int j = 0; j < allkmers[i].size(); ++j)
		{
			h_kmers[tmp++] = allkmers[i][j].getArray()[0];
			h_kmers[tmp++] = allkmers[i][j].getArray()[1];
			h_kmers[tmp++] = allkmers[i][j].getArray()[2];
			h_kmers[tmp++] = allkmers[i][j].getArray()[3];
		}
	}

	cudaProfilerStart();
	auto t1 = Clock::now();
	uint64_t *d_kmers = NULL;
	cudaMalloc((void **)&d_kmers, sizeof(*d_kmers) * totkmers * N_LONGS);
	cudaMemcpy(d_kmers, h_kmers, sizeof(*d_kmers) * totkmers * N_LONGS,
			   cudaMemcpyHostToDevice);

	// Bloom filter construction
	typedef nvbio::bloom_filter<5, RSHash<uint64_t *>,
								ElfHash<uint64_t *>, uint32_t *>
		bloom_filter_type;
	uint64_t nfilter_elems = totkmers / 4;
	uint32_t *d_filter_storage = NULL;
	cudaMalloc((void **)&d_filter_storage,
			   nfilter_elems * sizeof(*d_filter_storage));
	cudaMemset(d_filter_storage, 0,
			   nfilter_elems * sizeof(*d_filter_storage));
	bloom_filter_type d_filter(nfilter_elems * 32, d_filter_storage);
	cout << "number of bits " << nfilter_elems * 32
		 << " " << (nfilter_elems * 32) / ((1 << 20) * 8) << " mb " << endl;

	// for counting
	uint8_t		 *d_kmer_pass = NULL;
	uint64_t	**d_kmer_ptrs = NULL;
	cudaMalloc((void **)&d_kmer_pass,
			   totkmers * sizeof(*d_kmer_pass));
	cudaMemset(d_kmer_pass, 0,
			   totkmers * sizeof(*d_kmer_pass));
	cudaMalloc((void **)&d_kmer_ptrs,
			   totkmers * sizeof(*d_kmer_ptrs));

	int nblocks = (totkmers + 1023) / 1024;
	populate_kernel<<<nblocks,1024>>>(totkmers, d_kmers, d_filter,
									  d_kmer_pass, d_kmer_ptrs);
	cudaDeviceSynchronize();

	cudaProfilerStop();

	// uint32_t *h_filter_storage = (uint32_t *)
	// 	malloc(sizeof(*h_filter_storage) * nfilter_elems);	
	// cudaMemcpy(h_filter_storage, d_filter_storage,
	// 		   nfilter_elems * sizeof(*h_filter_storage),
	// 		   cudaMemcpyDeviceToHost);
	// bloom_filter_type h_filter(nfilter_elems * 32, h_filter_storage);

	// for (int i = 0; i < 2; ++i)
	// {
	// 	int kmerid = rand() % totkmers;
	// 	{
	// 		uint64_t *kmerptr = &(h_kmers[kmerid * N_LONGS]);
	// 		printf("%llu %llu %llu %llu\n",
	// 			   kmerptr[0], kmerptr[1], kmerptr[2], kmerptr[3]);
	// 		size_t i,j,l;
	// 		char *sx = (char *) malloc(1024);
	// 		char *s = sx;
	// 		memset(s, '\0', 1024);
	// 		printf("kmerid %d %d\n", kmerid, Kmer::k);
	// 		for (i = 0; i < Kmer::k; i++)
	// 		{
	// 			j = i % 32;
	// 			l = i / 32;

	// 			switch(((kmerptr[l]) >> (2*(31-j)) ) & 0x03)
	// 			{
	// 			case 0x00: *s = 'A'; ++s; break;
	// 			case 0x01: *s = 'C'; ++s; break;
	// 			case 0x02: *s = 'G'; ++s; break;
	// 			case 0x03: *s = 'T'; ++s; break;
	// 			}
	// 		}

	// 		bool b = h_filter.has(&(h_kmers[kmerid]));
	// 		printf("kmer %s in filter? %d\n", sx, b);
	// 	}
	// }
	
	printf("Sorting to remove filtered kmers...\n");
	thrust::device_ptr<uint8_t> d_kmer_pass_tht(d_kmer_pass);
	thrust::device_ptr<uint64_t *> d_kmer_ptrs_tht(d_kmer_ptrs);
	thrust::sort_by_key(d_kmer_pass_tht, d_kmer_pass_tht + totkmers,
						d_kmer_ptrs_tht, thrust::greater<uint8_t>());

	thrust::device_vector<uint8_t>::iterator iter;
	iter = thrust::find(d_kmer_pass_tht, d_kmer_pass_tht + totkmers, 0);

	uint8_t *x = thrust::raw_pointer_cast(&(iter[0]));
	int count = x - d_kmer_pass;
	
	// printf("%p %p --- active: %d, elim: %d, tot: %d\n",
	// 	   x, d_kmer_pass, count, totkmers-count, totkmers);

	// sort
	uint64_t active_kmers = x - d_kmer_pass;
	thrust::sort(d_kmer_ptrs_tht, d_kmer_ptrs_tht + active_kmers,
				 compare_kmer());

	printf("reducing...\n");

	// reduce by key
	uint64_t **d_outkeys = NULL;
	uint32_t *d_outvals = NULL;
	cudaMalloc((void **)&d_outkeys, active_kmers * sizeof(*d_outkeys));
	cudaMalloc((void **)&d_outvals, active_kmers * sizeof(*d_outvals));
	thrust::device_ptr<uint64_t *> d_outkeys_tht(d_outkeys);
	thrust::device_ptr<uint32_t> d_outvals_tht(d_outvals);

	typedef thrust::device_vector<uint64_t *>::iterator uint64PtrIter;
	typedef thrust::device_vector<uint32_t>::iterator uint32Iter;
	typedef thrust::pair<uint64PtrIter, uint32Iter> pairIter;
	pairIter var = 
		thrust::reduce_by_key(d_kmer_ptrs_tht, d_kmer_ptrs_tht + active_kmers,
							  d_kmer_pass_tht,
							  d_outkeys_tht, d_outvals_tht,
							  equal_kmer(),
							  thrust::plus<uint32_t>());
	uint64PtrIter	iter1 = thrust::get<0>(var);
	uint64_t		distinct_kmers =
		thrust::raw_pointer_cast(&(iter1[0])) - d_outkeys;

	cout << "distinct kmers: " << distinct_kmers << "\n";

	uint64_t *tmp2 = d_outkeys_tht[17];
	uint32_t x2 = d_outvals_tht[17];
	
	// printf("%p %p %d\n", d_kmers, tmp2, x2);

    uint64_t **h_outkeys = (uint64_t **)
		malloc(sizeof(*h_outkeys) * distinct_kmers);
	uint32_t *h_outvals = (uint32_t *)
		malloc(sizeof(*h_outvals) * distinct_kmers);

	cudaMemcpy(h_outkeys, d_outkeys,
			   sizeof(*h_outkeys) * distinct_kmers,
			   cudaMemcpyDeviceToHost);
	cudaMemcpy(h_outvals, d_outvals,
			   sizeof(*h_outvals) * distinct_kmers,
			   cudaMemcpyDeviceToHost);
	uint64_t *base_kmer_ptr = d_kmers;


	auto t2 = Clock::now();
	double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
	duration = duration / 1e6;
	printf("bloom filter insert/query took %.2f milliseconds on GPU for %d kmers\n",
		   duration, totkmers);


	// printf("%p %p %d\n", h_outkeys[0], base_kmer_ptr,
	// 	   (h_outkeys[0] - base_kmer_ptr)/N_LONGS);
	// for (uint64_t i = 0; i < distinct_kmers; ++i)
	// 	assert ((h_outkeys[i] - base_kmer_ptr)/N_LONGS < totkmers);

	// {
	// 	uint64_t kmerid = (h_outkeys[77777] - base_kmer_ptr)/N_LONGS;
	// 	uint64_t *kmerptr = &(h_kmers[kmerid * N_LONGS]);
	// 	printf("%llu %llu %llu %llu\n",
	// 		   kmerptr[0], kmerptr[1], kmerptr[2], kmerptr[3]);
	// 	size_t i,j,l;
	// 	char *sx = (char *) malloc(1024);
	// 	char *s = sx;
	// 	memset(s, '\0', 1024);
	// 	printf("kmerid %d %d\n", kmerid, Kmer::k);
	// 	for (i = 0; i < Kmer::k; i++)
	// 	{
	// 		j = i % 32;
    // 		l = i / 32;
	// 		printf("%d\n", i);

	// 		switch(((kmerptr[l]) >> (2*(31-j)) ) & 0x03)
	// 		{
	// 		case 0x00: *s = 'A'; ++s; break;
	// 		case 0x01: *s = 'C'; ++s; break;
	// 		case 0x02: *s = 'G'; ++s; break;
	// 		case 0x03: *s = 'T'; ++s; break;
	// 		}
	// 	}

	// 	printf("kmer %s occurs %d\n", sx, h_outvals[77777]);
	// }

	

	// ptrdiff_t *d_outids = NULL;
	// cudaMalloc((void **)&d_outids, distinct_kmers * sizeof(*d_outids));
	// thrust::device_ptr<ptrdiff_t> d_outids_tht(d_outids);
	// thrust::transform(d_outkeys_tht, d_outkeys_tht + distinct_kmers,
	// 				  thrust::constant_iterator<uint64_t *>(d_kmers),
	// 				  d_outids_tht,
	// 				  minus_ptr<uint64_t *>());

	// ptrdiff_t x3 = d_outids[0];
	// cout << x3 << "\n";
	
	
	
	// uint32_t count = h_ptr - d_kmer_pass;
	// cout << h_ptr << " " << d_kmer_pass
	// 	 <<  " # of passed kmers " << count << "\n";	


	// uint32_t *h_filter_storage = (uint32_t *)
	// 	malloc(sizeof(*h_filter_storage) * nfilter_elems);	
	// cudaMemcpy(h_filter_storage, d_filter_storage,
	// 		   nfilter_elems * sizeof(*h_filter_storage),
	// 		   cudaMemcpyDeviceToHost);
	// bloom_filter_type h_filter(nfilter_elems * 32, h_filter_storage);

	// uint64_t tmpx[4]; tmpx[0] = 2131321; tmpx[1] = 21312; tmpx[2] = 99; tmpx[3] = 123123;
	// std::cout << h_filter.has(&(h_kmers[6])) << "\n";
	// std::cout << h_filter.has(tmpx) << "\n";
	////////////////////////////////////////////////////////////////////////////
	

}
#endif
