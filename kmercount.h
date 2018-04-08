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

#include "libcuckoo/cuckoohash_map.hh"
#include "libbloom/bloom64.h"

#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include "kmercode/Buffer.h"
#include "kmercode/common.h"
#include "kmercode/fq_reader.h"
#include "kmercode/ParallelFASTQ.h"
#include "kmercode/rbounds.hpp"
#include "kmercode/hyperloglog.hpp"

typedef cuckoohash_map<Kmer, size_t> dictionary_t; // <k-mer && reverse-complement, #kmers>

struct filedata {

    char filename[MAX_FILE_PATH];
    size_t filesize;
};

/**
 * @brief GetFiles
 * @param filename
 * @return
 */
std::vector<filedata>  GetFiles(char *filename) {
    int64_t totalsize = 0;
    int numfiles = 0;
    std::vector<filedata> filesview;
    
    filedata fdata;
    ifstream allfiles(filename);
    if(!allfiles.is_open()) {
        cerr << "could not open " << filename << endl;
        exit(1);
    }
    allfiles.getline(fdata.filename,MAX_FILE_PATH);
    while(!allfiles.eof())
    {
        struct stat st;
        stat(fdata.filename, &st);
        fdata.filesize = st.st_size;
        
        filesview.push_back(fdata);
        cout << filesview.back().filename << " : " << filesview.back().filesize / (1024*1024) << " MB" << endl;
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
    double kdict = omp_get_wtime();
	
    ifstream filein(kmer_file);
    string line;
    size_t elem;
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
                kmerfromstr.set_kmer(kmerstr.c_str());

                auto updatecountjelly = [&elem](size_t &num) { num+=elem; };
                // If the number is already in the table, it will increment its count by the occurrence of the new element. 
                // Otherwise it will insert a new entry in the table with the corresponding k-mer occurrence.
                countsjelly.upsert(kmerfromstr.rep(), updatecountjelly, elem);      
            }
    } else std::cout << "Unable to open the input file\n";
    filein.close();
    //cout << "jellyfish file parsing took: " << omp_get_wtime()-kdict << "s" << endl;
    //
    // Reliable k-mer filter on countsjelly
    //
    int kmer_id = 0;

    auto ltj = countsjelly.lock_table(); // jellyfish counting
    for (const auto &it : ltj) 
        if (it.second >= lower && it.second <= upper)
        {
            countsreliable_jelly.insert(it.first,kmer_id);
            ++kmer_id;
        }
    ltj.unlock(); // unlock the table
    // Print some information about the table
    //cout << "Table size Jellyfish: " << countsjelly.size() << std::endl;
    //cout << "Entries within reliable range Jellyfish: " << countsreliable_jelly.size() << std::endl;    
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
void DeNovoCount(vector<filedata> & allfiles, dictionary_t & countsreliable_denovo, int lower, int upper, int kmer_len, size_t upperlimit /* memory limit */)
{
    vector < vector<Kmer> > allkmers(MAXTHREADS);
    vector < HyperLogLog > hlls(MAXTHREADS, HyperLogLog(12));	// std::vector fill constructor    
    
    double denovocount = omp_get_wtime();
    dictionary_t countsdenovo;
    //cout << "\nRunning with up to " << MAXTHREADS << " threads" << endl;
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
		size_t tlreads = 0;	// thread local reads

        	size_t fillstatus = 1;
        	while(fillstatus) 
		{ 
            		fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
            		size_t nreads = seqs.size();
            
       	    		//#pragma omp parallel for
			for(int i=0; i<nreads; i++) 
			{
                		// remember that the last valid position is length()-1
                		int len = seqs[i].length();
           
                		for(int j=0; j<=len-kmer_len; j++)  
                		{
                    			std::string kmerstrfromfastq = seqs[i].substr(j, kmer_len);
                   			Kmer mykmer(kmerstrfromfastq.c_str());
                    			Kmer lexsmall = mykmer.rep();  
		    			allkmers[MYTHREAD].push_back(lexsmall);
		   			hlls[MYTHREAD].add((const char*) lexsmall.getBytes(), lexsmall.getNumBytes());   
				}
			} // for(int i=0; i<nreads; i++)
			tlreads += nreads;
            	} //while(fillstatus) 
		delete pfq;

		#pragma omp critical
		totreads += tlreads;
	}
	//cout << "There were " << totreads << " reads" << endl;
    }  

    // HLL reduction (serial for now)
    for (int i=1;i< MAXTHREADS; i++) 
    {
	    std::transform(hlls[0].M.begin(), hlls[0].M.end(), hlls[i].M.begin(), hlls[0].M.begin(), [](uint8_t c1, uint8_t c2) -> uint8_t{ return std::max(c1, c2); });
    }
    double cardinality = hlls[0].estimate();
    //cout << "Cardinality estimate is " << cardinality << endl;


    unsigned int random_seed = 0xA57EC3B2;
    const double desired_probability_of_false_positive = 0.05;
    struct bloom * bm = (struct bloom*) malloc(sizeof(struct bloom));
    bloom_init64(bm, cardinality * 1.1, desired_probability_of_false_positive);
    //cout << "Table size is: " << bm->bits << " bits, " << ((double)bm->bits)/8/1024/1024 << " MB" << endl;
    //cout << "Optimal number of hash functions is : " << bm->hashes << endl;

    #pragma omp parallel
    {	    
	for(auto v:allkmers[MYTHREAD])
	{
		// ABAB: race condition here
		// thread-1 can start checking and find a zero bit in kth hash for item t, 
		// but then thread-2 takes over and finds a zero bit in (k+i)th hash for the same item t
		// both threads will think bloom filter never saw t yet, and will not insert into the 
		// hash table; potentially causing the item (k-mer) t to be lost. 
		bool inBloom = (bool) bloom_check_add(bm, v.getBytes(), v.getNumBytes(),1);

		if(inBloom) countsdenovo.insert(v, 0);
	}
    }


    // in this pass, only use entries that already are in the hash table
    auto updatecount = [](size_t &num) { ++num; };
    #pragma omp parallel
    {	    
	for(auto v:allkmers[MYTHREAD])
	{
		// does nothing if the entry doesn't exist in the table
		countsdenovo.update_fn(v,updatecount);
	}
    }

    cout << "\ndenovo counting took: " << omp_get_wtime()-denovocount << "s" << endl;


    //
    // Reliable k-mer filter on countsdenovo
    //
    size_t kmer_id_denovo = 0;
    
    auto lt = countsdenovo.lock_table(); // our counting
    for (const auto &it : lt) 
        if (it.second >= lower && it.second <= upper)
        {
            countsreliable_denovo.insert(it.first,kmer_id_denovo);
            ++kmer_id_denovo;
        }
    lt.unlock(); // unlock the table
    // Print some information about the table
    //cout << "Table size: " << countsdenovo.size() << std::endl;
    cout << "Entries within reliable range: " << countsreliable_denovo.size() << std::endl;    
    //cout << "Bucket count: " << countsdenovo.bucket_count() << std::endl;
    //cout << "Load factor: " << countsdenovo.load_factor() << std::endl;
    countsdenovo.clear(); // free  
 
}
#endif
