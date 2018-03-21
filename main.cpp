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
//#define JELLYFISH

using namespace std;

#ifdef _ALLKMER
struct spmatType_ {

    int count = 0;   /* number of shared k-mers */
    vector<std::pair<int,int>> vpos; /* wanna keep all the positions */
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
typedef cuckoohash_map<Kmer, int> dictionary_t; // <k-mer && reverse-complement, #kmers>
typedef std::vector<Kmer> Kmers;

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
#ifdef JELLYFISH
                if(thisOpt->argument == NULL)
                {
                    cout << "BELLA execution terminated: -f requires an argument" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
#endif
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
                cout << " -f : k-mer list from Jellyfish (required if Jellyfish k-mer counting is used)" << endl; // the reliable k-mers are selected by bella
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

#ifdef JELLYFISH
    if(kmer_file == NULL || all_inputs_fofn == NULL || out_file == NULL || depth == 0)
    {
        cout << "BELLA execution terminated: missing arguments" << endl;
        cout << "Run with -h to print out the command line options\n" << endl;
        return 0;
    }
#else
    if(all_inputs_fofn == NULL || out_file == NULL || depth == 0)
    {
        cout << "BELLA execution terminated: missing arguments" << endl;
        cout << "Run with -h to print out the command line options\n" << endl;
        return 0;
    }
#endif

    free(optList);
    free(thisOpt);

    //
    // Declarations 
    //
#ifdef JELLYFISH
    ifstream filein(kmer_file);
#endif
    vector<filedata> allfiles = GetFiles(all_inputs_fofn);
    FILE *fastafile;
    int elem;
    int lower = 2; // reliable range lower bound (fixed)
    int upper;     // reliable range upper bound
    char *buffer;
    string kmerstr;
    string line;
    Kmer::set_k(kmer_len);
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
#ifdef JELLYFISH
    cout << "k-mer counting: Jellyfish" << endl;
    cout << "Input k-mer file: " << kmer_file << endl;
#endif
    cout << "k-mer counting: BELLA" << endl;
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
    double all = omp_get_wtime();
    upper = rbounds(depth,erate,kmer_len);
    cout << "Reliable lower bound: " << lower << endl;
    cout << "Reliable upper bound: " << upper << "\n" << endl;

    //
    // Reliable k-mer file parsing and k-mer dictionary creation
    // NOTE: this will be replaced by our k-mer counting
    //
    //
#ifdef JELLYFISH
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
    dictionary_t countsreliable_jelly; 
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
#else
    double denovocount = omp_get_wtime();
    dictionary_t countsdenovo;

    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);

        size_t fillstatus = 1;
        while(fillstatus) { 
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
    cout << "\ndenovo counting took: " << omp_get_wtime()-denovocount << "s" << endl;
    //
    // Reliable k-mer filter on countsdenovo
    //
    dictionary_t countsreliable_denovo; 
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
#endif   
    //
    // Fastq(s) parsing
    //
    double parsefastq = omp_get_wtime();
    int read_id = 0; // read_id needs to be global (not just per file)

#ifdef _OPENMP
    int numthreads = omp_get_max_threads();
    cout << "\nRunning with up to " << numthreads << " threads" << endl;
#else
    int numthreads = 1;
#endif       
        vector < vector<tuple<int,int,int>> > alloccurrences(numthreads);   
        vector < vector<tuple<int,int,int>> > alltranstuples(numthreads);   
        vector < readVector_ > allreads(numthreads);

        double time1 = omp_get_wtime();

        for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);

        size_t fillstatus = 1;
        while(fillstatus) { 
            fillstatus = pfq->fill_block(nametags, seqs, quals, upperlimit);
            int nreads = seqs.size();
       
                #pragma omp parallel for
                for(int i=0; i<nreads; i++) 
                {
                    // remember that the last valid position is length()-1
                    int len = seqs[i].length();
                    
                    readType_ temp;
                    temp.nametag = nametags[i];
                    temp.seq = seqs[i]; // save reads for seeded alignment
                    temp.readid = read_id+i;
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
        
                        int idx; // kmer_id
            #ifdef JELLYFISH
                        auto found = countsreliable_jelly.find(lexsmall,idx);
            #else
                        auto found = countsreliable_denovo.find(lexsmall,idx);
            #endif
                        if(found)
                        {
                            alloccurrences[tid].emplace_back(std::make_tuple(read_id+i,idx,j)); // vector<tuple<read_id,kmer_id,kmerpos>>
                            alltranstuples[tid].emplace_back(std::make_tuple(idx,read_id+i,j)); // transtuples.push_back(col_id,row_id,kmerpos)
                        }
                    }
                }   // for(int i=0; i<nreads; i++)
                    //cout << "total number of reads processed so far is " << read_id << endl;
            read_id += nreads;
        } //while(fillstatus) 
    delete pfq;
    } // for all files
    cout << "global vector filling took: " << omp_get_wtime()-time1 << "s" << endl;

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
    cout << "\nTotal number of reads: "<< read_id << endl;
    cout << "fastq(s) parsing fastq took: " << omp_get_wtime()-parsefastq << "s" << endl;
    double matcreat = omp_get_wtime();

    //
    // Sparse matrices construction
    //

#ifdef JELLYFISH
    int nkmer = countsreliable_jelly.size();
#else
    int nkmer = countsreliable_denovo.size();
#endif

    CSC<int, int> spmat(occurrences, read_id, nkmer, 
                            [] (int & p1, int & p2) 
                            {   // assume no errors in MergeDuplicates
                                // keep just the first position of that k-mer in that read
                                return p1;
                            });
    //std::cout << "spmat created with " << spmat.nnz << " nonzeros" << endl;
    std::vector<tuple<int,int,int>>().swap(occurrences);    // remove memory of occurences

    CSC<int, int> transpmat(transtuples, nkmer, read_id, 
                            [] (int & p1, int & p2) 
                            {  return p1;
                            });
    //std::cout << "transpose(spmat) created" << endl;
    std::vector<tuple<int,int,int>>().swap(transtuples); // remove memory of transtuples

    spmat.Sorted();
    transpmat.Sorted();
    cout << "spmat and spmat^T creation took: " << omp_get_wtime()-matcreat << "s" << endl;
    
    //
    // Overlap detection (sparse matrix multiplication) and seed-and-extend alignment
    //

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
    cout << "total running time: " << omp_get_wtime()-all << "s\n" << endl;
    return 0;
} 