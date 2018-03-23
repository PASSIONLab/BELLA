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

#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include "kmercode/Buffer.h"
#include "kmercode/common.h"
#include "kmercode/fq_reader.h"
#include "kmercode/ParallelFASTQ.h"
#include "kmercode/rbounds.hpp"
#include "kmercode/hyperloglog.hpp"

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
                kmerfromstr.set_kmer(kmerstr.c_str());

                auto updatecountjelly = [&elem](int &num) { num+=elem; };
                // If the number is already in the table, it will increment its count by the occurrence of the new element. 
                // Otherwise it will insert a new entry in the table with the corresponding k-mer occurrence.
                countsjelly.upsert(kmerfromstr.rep(), updatecountjelly, elem);      
            }
    } else std::cout << "Unable to open the input file\n";
    filein.close();
    cout << "jellyfish file parsing took: " << omp_get_wtime()-kdict << "s" << endl;
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
    cout << "Table size Jellyfish: " << countsjelly.size() << std::endl;
    cout << "Entries within reliable range Jellyfish: " << countsreliable_jelly.size() << std::endl;    
    cout << "Bucket count Jellyfish: " << countsjelly.bucket_count() << std::endl;
    cout << "Load factor Jellyfish: " << countsjelly.load_factor() << std::endl;
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
    vector<string> seqs;
    vector<string> quals;
    vector<string> nametags;
    vector<Kmer> allkmers;
    HyperLogLog hll(12);
    
    
    double denovocount = omp_get_wtime();
    dictionary_t countsdenovo;

    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) 
    {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);

        size_t fillstatus = 1;
        while(fillstatus) 
	{ 
            fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
            int nreads = seqs.size();
	    auto updatecount = [](int &num) { ++num; };
            
       	    #pragma omp parallel for
            for(int i=0; i<nreads; i++) 
            {
                // remember that the last valid position is length()-1
                int len = seqs[i].length();

#ifdef _OPENMP
		int tid = omp_get_thread_num();
#else
		int tid = 0;
#endif
           
                for(int j=0; j<=len-kmer_len; j++)  
                {
                    std::string kmerstrfromfastq = seqs[i].substr(j, kmer_len);
                    Kmer mykmer(kmerstrfromfastq.c_str());
                    // remember to use only ::rep() when building kmerdict as well
                    Kmer lexsmall = mykmer.rep();  
		    //allkmers.push_back(lexsmall);

		    hll.add((const char*) lexsmall.getBytes(), lexsmall.getNumBytes());   //ABAB: trivial to make this multithreaded, but what about bloom filter?
		    
		    // If the number is already in the table, it will increment its count by one. 
		    // Otherwise it will insert a new entry in the table with count one.
		    countsdenovo.upsert(lexsmall, updatecount, 1);
		}
            } // for(int i=0; i<nreads; i++)
	} //while(fillstatus) 
	delete pfq;
    }  
    cout << "\ndenovo counting took: " << omp_get_wtime()-denovocount << "s" << endl;
    //
    // Reliable k-mer filter on countsdenovo
    //
    int kmer_id_denovo = 0;
    
    auto lt = countsdenovo.lock_table(); // our counting
    for (const auto &it : lt) 
        if (it.second >= lower && it.second <= upper)
        {
            countsreliable_denovo.insert(it.first,kmer_id_denovo);
            ++kmer_id_denovo;
        }
    lt.unlock(); // unlock the table
    // Print some information about the table
    cout << "Table size: " << countsdenovo.size() << std::endl;
    cout << "Entries within reliable range: " << countsreliable_denovo.size() << std::endl;    
    cout << "Bucket count: " << countsdenovo.bucket_count() << std::endl;
    cout << "Load factor: " << countsdenovo.load_factor() << std::endl;
    countsdenovo.clear(); // free  
 
    double cardinality = hll.estimate();
    cout << "Cardinality estimate is " << cardinality << endl;
}
#endif
