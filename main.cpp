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

#include "mtspgemm2017/utility.h"
#include "mtspgemm2017/CSC.h"
#include "mtspgemm2017/CSR.h"
#include "mtspgemm2017/global.h"
#include "mtspgemm2017/IO.h"
#include "mtspgemm2017/multiply.h"


#define LSIZE 16000
#define ITERS 10
#define DENOVOCOUNT
//#define DEPTH 30
//#define _ALLKMER


using namespace std;

#ifdef _ALLKMER
struct spmatType_ {

    int count = 0;   /* number of shared k-mers */
    std::vector<std::pair<int,int>> vpos; /* wanna keep all the positions */
};
#else
struct spmatType_ {

    int count = 0;   /* number of shared k-mers */
    int pos[4] = {0};  /* pos1i, pos1j, pos2i, pos2j */
};
#endif

struct filedata {

    char filename[MAX_FILE_PATH];
    size_t filesize;
};

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

typedef shared_ptr<spmatType_> spmatPtr_; // pointer to spmatType_ datastruct
// typedef std::unordered_map<Kmer, int> dictionary; // <k-mer && reverse-complement, #kmers>
typedef cuckoohash_map<Kmer, int> dictionary; // <k-mer && reverse-complement, #kmers>
typedef std::vector<Kmer> Kmers;

// Function to create the dictionary
// assumption: kmervect has unique entries
void dictionaryCreation(dictionary &kmerdict, Kmers &kmervect)
{
    kmerdict.reserve(kmervect.size());	// only makes sense for std::unordered_map
    for(int i = 0; i<kmervect.size(); i++)
    {
        // kmerdict.insert(make_pair(kmervect[i].rep(), i));
        kmerdict.insert(kmervect[i].rep(), i);
	
    }
}

int main (int argc, char *argv[]) {

    //
    // Program name and purpose
    //
    cout << "\nBELLA - Long Read Aligner for De Novo Genome Assembly\n" << endl;
    //
    // Setup the input files
    //

    option_t *optList, *thisOpt;
    // Get list of command line options and their arguments 
    optList = NULL;
    optList = GetOptList(argc, argv, (char*)"f:i:o:d:hk:a:ze:p:");
   
    bool skip_algnmnt_krnl = false;
    char *kmer_file = NULL; // Reliable k-mer file from Jellyfish
    char *all_inputs_fofn = NULL;   // List of fastq(s)
    char *out_file = NULL; // output filename
    int kmer_len = 17;  // default k-mer length
    int algnmnt_thr = 50;   // default alignment score threshold
    int algnmnt_drop = 3;   // default alignment x-drop factor
    double erate = 0.15; // default error rate
    int depth = 0; // depth/coverage required

    if(optList == NULL)
    {
        cout << "BELLA execution terminated: not enough parameters or invalid option" << endl;
        cout << "Run with -h to print out the command line options\n" << endl;
        return 0;
    }

    while (optList!=NULL) {
        thisOpt = optList;
        optList = optList->next;
        switch (thisOpt->option) {
            case 'f': {
                if(thisOpt->argument == NULL)
                {
                    cout << "BELLA execution terminated: -f requires an argument" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                kmer_file = strdup(thisOpt->argument);
                break;
            }
            case 'i': {
                if(thisOpt->argument == NULL)
                {
                    cout << "BELLA execution terminated: -i requires an argument" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                all_inputs_fofn = strdup(thisOpt->argument);
                break;
            }
            case 'o': {
                if(thisOpt->argument == NULL)
                {
                    cout << "BELLA execution terminated: -o requires an argument" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                char* line1 = strdup(thisOpt->argument);
                char* line2 = strdup(".bla");
                size_t len1 = strlen(line1);
                size_t len2 = strlen(line2);

                out_file = (char*)malloc(len1 + len2 + 1);
                if (!out_file) abort();

                memcpy(out_file, line1, len1);
                memcpy(out_file + len1, line2, len2);
                out_file[len1 + len2] = '\0';

                delete line1;
                delete line2;

                break;
            }
            case 'd': {
                if(thisOpt->argument == NULL)
                {
                    cout << "BELLA execution terminated: -d requires an argument" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                depth = atoi(thisOpt->argument);  
                break;
            }
            case 'z': skip_algnmnt_krnl = true; break; // TO DO: add skip alignment
            case 'k': {
                kmer_len = atoi(thisOpt->argument);
                break;
            }
            case 'e': {
                erate = strtod(thisOpt->argument, NULL);
                break;
            }
            case 'a': {
                algnmnt_thr = atoi(thisOpt->argument);
                break;
            }
            case 'p': {
                algnmnt_drop = atoi(thisOpt->argument);
                break;
            }
            case 'h': {
                cout << "Usage:\n" << endl;
                cout << " -f : k-mer list from Jellyfish (required)" << endl; // the reliable k-mers are selected by bella
                cout << " -i : list of fastq(s) (required)" << endl;
                cout << " -o : output filename (required)" << endl;
                cout << " -d : depth/coverage (required)" << endl; // TO DO: add depth estimation
                cout << " -k : k-mer length [17]" << endl;
                cout << " -a : alignment score threshold [50]" << endl;
                cout << " -p : alignment x-drop factor [3]" << endl;
                cout << " -e : error rate [0.15]" << endl;
                cout << " -z : skip the alignment [false]\n" << endl;
                
                FreeOptList(thisOpt); // Done with this list, free it
                return 0;
            }
        }
    }

    if(kmer_file == NULL || all_inputs_fofn == NULL || out_file == NULL || depth == 0)
    {
        cout << "BELLA execution terminated: missing arguments" << endl;
        cout << "Run with -h to print out the command line options\n" << endl;
        return 0;
    }
  
    free(optList);
    free(thisOpt);

    //
    // Declarations 
    //
    ifstream filein(kmer_file);
    vector<filedata> allfiles = GetFiles(all_inputs_fofn);
    FILE *fastafile;
    int elem;
    int lower = 2; // reliable range lower bound (fixed)
    int upper;     // reliable range upper bound
    char *buffer;
    string kmerstr;
    string line;
    Kmer::set_k(kmer_len);
    dictionary kmerdict;
    size_t upperlimit = 10000000; // in bytes
    Kmer kmerfromstr;
    Kmers kmervect;
    vector<string> seqs;
    vector<string> quals;
    vector<string> nametags;
    readVector_ reads;
    Kmers kmersfromreads;
    vector<tuple<int,int,int>> occurrences;
    vector<tuple<int,int,int>> transtuples;
    
    // 
    // File and setting used
    //
    cout << "Input k-mer file: " << kmer_file << endl;
    cout << "Output file: " << out_file << endl;
    cout << "k-mer length: " << kmer_len << endl;
    if(skip_algnmnt_krnl)
        cout << "Compute alignment: false" << endl;
    else cout << "Compute alignment: true" << endl;
    cout << "Alignment x-drop factor: " << algnmnt_drop << endl;
    cout << "Alignment score threshold: " << algnmnt_thr << endl;
    cout << "Depth: " << depth << "X" << endl;

    //
    // Reliable bounds computation
    //
    upper = rbounds(depth,erate,kmer_len);

    cout << "Reliable lower bound: " << lower << endl;
    cout << "Reliable upper bound: " << upper << "\n" << endl;

    //
    // Reliable k-mer file parsing and k-mer dictionary creation
    // NOTE: this will be replaced by our k-mer counting
    //
    //

#ifdef DENOVOCOUNT
    double denovocount = omp_get_wtime();
    dictionary countsdenovo;

    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);

        size_t fillstatus = 1;
        while(fillstatus) { 
            double time1 = omp_get_wtime();
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

		    // If the number is already in the table, it will increment its count by one. 
		    // Otherwise it will insert a new entry in the table with count one.
		    countsdenovo.upsert(lexsmall, updatecount, 1);

		}
            } // for(int i=0; i<nreads; i++)
	} //while(fillstatus) 
    	delete pfq;
    }  
    cout << "denovo counting took: " << omp_get_wtime()-denovocount << "s\n" << endl;
    // Print some information about the table
    
    int64_t keep = 0;
    vector < pair<Kmer, int> > tokeep;
    vector < pair<Kmer, int> > verify;	// Aydin: remove after check
    vector < pair<Kmer, int> > errors;  // Aydin: remove after check
    
    
    auto lt = countsdenovo.lock_table();
    for (const auto &it : lt) {
      if (it.second >= lower && it.second <= upper)
	      tokeep.push_back(make_pair(it.first, it.second));      	
    }
    sort(tokeep.begin(), tokeep.end());

    cout << "Table size: " << countsdenovo.size() << std::endl;
    cout << "Entries within reliable range: " << keep << std::endl;    
    cout << "Bucket count: " << countsdenovo.bucket_count() << std::endl;
    cout << "Load factor: " << countsdenovo.load_factor() << std::endl;
    
#endif

    double all = omp_get_wtime();
    double kdict = omp_get_wtime();
    if(filein.is_open()) {
            while(getline(filein, line)) {
                if(line.length() == 0)
                    break;

                string substring = line.substr(1);
                elem = stoi(substring);
                getline(filein, kmerstr);   
                kmerfromstr.set_kmer(kmerstr.c_str());
                if(elem >= lower && elem <= upper) // keep just the kmers within the reliable bounds
		{
                    kmervect.push_back(kmerfromstr); 

                    Kmer mykmer(kmerstr.c_str());	// Aydin: remove after check	        
	    	    verify.push_back(make_pair (mykmer.rep(), elem) ); // Aydin: remove after check
		}		    
            }
    } else std::cout << "Unable to open the input file\n";
    filein.close();

    /* <begin> Aydin: remove after check */
    sort(verify.begin(), verify.end()); 
    set_difference( tokeep.begin(), tokeep.end(), verify.begin(), verify.end(), back_inserter( errors ) );     for(const auto & it : errors)
    {
	    cout << it.first.toString() << " " << it.second << endl;
    }
    /* <end> Aydin: remove after check */


    dictionaryCreation(kmerdict, kmervect);
    cout << "Reliable k-mer: " << kmerdict.size() << endl;
    cout << "k-mer dictionary creation took: " << omp_get_wtime()-kdict << "s\n" << endl;
    
    //
    // Fastq(s) parsing
    //
    double parsefastq = omp_get_wtime();
    int read_id = 0; // read_id needs to be global (not just per file)


#ifdef _OPENMP
	int numthreads = omp_get_max_threads();
	cout << "Running with up to " << numthreads << " threads" << endl;
#else
	int numthreads = 1;
#endif


       
        vector < vector<tuple<int,int,int>> > alloccurrences(numthreads);	
        vector < vector<tuple<int,int,int>> > alltranstuples(numthreads);	
	vector < readVector_ > allreads(numthreads);	

    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);

        size_t fillstatus = 1;
        while(fillstatus) { 
            double time1 = omp_get_wtime();
            fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
            int nreads = seqs.size();
            
	    #pragma omp parallel for
            for(int i=0; i<nreads; i++) 
            {
                // remember that the last valid position is length()-1
                int len = seqs[i].length();

        	readType_ temp;
                temp.nametag = nametags[i];
                temp.seq = seqs[i];
                temp.readid = read_id+i;
                // save reads for seeded alignment
#ifdef _OPENMP
		int tid = omp_get_thread_num();
#else
		int tid = 0;
#endif
                allreads[tid].push_back(temp);
                
                for(int j=0; j<=len-kmer_len; j++)  
                {
                    std::string kmerstrfromfastq = seqs[i].substr(j, kmer_len);
                    Kmer mykmer(kmerstrfromfastq.c_str());
                    // remember to use only ::rep() when building kmerdict as well
                    Kmer lexsmall = mykmer.rep();      

		    int out;
                    // auto found = kmerdict.find(lexsmall);
		    auto found = kmerdict.find(lexsmall, out);
		    if(found)
                    {
                        alloccurrences[tid].emplace_back(std::make_tuple(read_id+i,out, j)); // vector<tuple<read_id,kmer_id,kmerpos>>
                        alltranstuples[tid].emplace_back(std::make_tuple(out, read_id+i, j)); // transtuples.push_back(col_id, row_id, kmerpos)
                    }
                }
            } // for(int i=0; i<nreads; i++)
	    read_id += nreads;
            // cout << "processed reads in " << omp_get_wtime()-time2 << " seconds "<< endl; 
            // cout << "total number of reads processed so far is " << read_id << endl;

        } //while(fillstatus) 
    delete pfq;
    } // for all files

    size_t readcount = 0;
	size_t tuplecount = 0;
	for(int t=0; t<numthreads; ++t)
	{
		readcount += allreads[t].size();
		tuplecount += alloccurrences[t].size();
	}
	reads.resize(readcount);
	occurrences.resize(tuplecount);
	transtuples.resize(tuplecount);

	size_t readssofar = 0;
	size_t tuplesofar = 0;
	for(int t=0; t<numthreads; ++t)
	{
		copy(allreads[t].begin(), allreads[t].end(), reads.begin()+readssofar);
		readssofar += allreads[t].size();

		copy(alloccurrences[t].begin(), alloccurrences[t].end(), occurrences.begin() + tuplesofar);
		copy(alltranstuples[t].begin(), alltranstuples[t].end(), transtuples.begin() + tuplesofar);
		tuplesofar += alloccurrences[t].size();
	}

    std::sort(reads.begin(), reads.end());   // bool operator in global.h: sort by readid

    std::vector<string>().swap(seqs);        // free memory of seqs  
    std::vector<string>().swap(quals);       // free memory of quals
    cout << "Total number of reads: "<< read_id << endl;
    cout << "Parsing fastq took: " << omp_get_wtime()-parsefastq << "s\n" << endl;
    double matcreat = omp_get_wtime();

    //
    // Sparse matrices construction
    //

    #ifdef _ALLKMER
    CSC<int, int> spmat(occurrences, read_id, kmervect.size(), 
                            [] (int & p1, int & p2) 
                            {   // assume no errors in MergeDuplicates
                                // keep just the first position of that k-mer in that read
                                return p1;
                            });
    //std::cout << "spmat created with " << spmat.nnz << " nonzeros" << endl;
    std::vector<tuple<int,int,int>>().swap(occurrences);    // remove memory of occurences

    CSC<int, int> transpmat(transtuples, kmervect.size(), read_id, 
                            [] (int & p1, int & p2) 
                            {  return p1;
                            });
    //std::cout << "transpose(spmat) created" << endl;
    std::vector<tuple<int,int,int>>().swap(transtuples); // remove memory of transtuples
    #else
    CSC<int, int> spmat(occurrences, read_id, kmervect.size(), 
                            [] (int & p1, int & p2) 
                            {   // assume no errors in MergeDuplicates
                                // keep just the first position of that k-mer in that read
                                return p1;
                            });
    //std::cout << "spmat created with " << spmat.nnz << " nonzeros" << endl;
    std::vector<tuple<int,int,int>>().swap(occurrences);    // remove memory of occurences

    CSC<int, int> transpmat(transtuples, kmervect.size(), read_id, 
                            [] (int & p1, int & p2) 
                            {  return p1;
                            });
    //std::cout << "transpose(spmat) created" << endl;
    std::vector<tuple<int,int,int>>().swap(transtuples); // remove memory of transtuples
    #endif

    spmat.Sorted();
    transpmat.Sorted();
    cout << "Spmat and Spmat^T creation took: " << omp_get_wtime()-matcreat << "s\n" << endl;
    
    //
    // Overlap detection (sparse matrix multiplication) and seed-and-extend alignment
    //

    #ifdef _ALLKMER
    spmatPtr_ getvaluetype(make_shared<spmatType_>());
    HeapSpGEMM(spmat, transpmat, 
            [] (int & pi, int & pj) // n-th k-mer positions on read i and on read j 
            {   spmatPtr_ value(make_shared<spmatType_>());
                value->count = 1;
                value->vpos.push_back(make_pair(pi,pj));
                return value;
            }, 
            [] (spmatPtr_ & m1, spmatPtr_ & m2)
            {   m2->count = m1->count+m2->count;
                // insert control on independent k-mer
                m2->vpos.insert(m2->vpos.end(), m1->vpos.begin(), m1->vpos.end());
                return m2;
            }, reads, getvaluetype, kmer_len, algnmnt_drop, algnmnt_thr, out_file, skip_algnmnt_krnl);
    #else
    spmatPtr_ getvaluetype(make_shared<spmatType_>());
    HeapSpGEMM(spmat, transpmat, 
            [] (int & pi, int & pj) // n-th k-mer positions on read i and on read j 
            {   spmatPtr_ value(make_shared<spmatType_>());
                value->count = 1;
                value->pos[0] = pi; // row
                value->pos[1] = pj; // col
                return value;
            }, 
            [] (spmatPtr_ & m1, spmatPtr_ & m2)
            {   m2->count = m1->count+m2->count;
                // m1->pos[0] = m1->pos[0]; // row 
                // m1->pos[1] = m1->pos[1]; // col
                m2->pos[2] = m1->pos[0]; // row 
                m2->pos[3] = m1->pos[1]; // col
                return m2;
            }, reads, getvaluetype, kmer_len, algnmnt_drop, algnmnt_thr, out_file, skip_algnmnt_krnl);
    #endif 
    cout << "Total running time: " << omp_get_wtime()-all << "s\n" << endl;
    return 0;
} 
