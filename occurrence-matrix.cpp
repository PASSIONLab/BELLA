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
#include <stack>
#include <cstring>
#include <string.h>
#include <cassert>
#include <ios>
#include <unordered_map>
#include <sys/stat.h>
#include <map>
#include <omp.h>
#include "kmercode/hash_funcs.h"
#include "kmercode/Kmer.hpp"
#include "kmercode/Buffer.h"
#include "kmercode/common.h"
#include "kmercode/fq_reader.h"
#include "kmercode/ParallelFASTQ.h"

#include "mtspgemm2017/utility.h"
#include "mtspgemm2017/CSC.h"
#include "mtspgemm2017/CSR.h"
#include "mtspgemm2017/IO.h"
#include "mtspgemm2017/multiply.h"

#define LSIZE 16000
#define KMER_LENGTH 15
#define ITERS 10

using namespace std;

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

typedef std::map<Kmer,size_t> dictionary; // <k-mer && reverse-complement, #kmers>
typedef std::vector<Kmer> Kmers;
typedef std::pair<size_t, vector<size_t> > cellspmat; // pair<kmer_id_j, vector<posix_in_read_i>>
typedef std::map<size_t, std::pair< std::vector<size_t>, std::vector<size_t> > > multcell; // map<kmer_id, vector<posix_in_read_i, posix_in_read_j>>

// Function to create the dictionary
// assumption: kmervect has unique entries
void dictionaryCreation(dictionary &kmerdict, Kmers &kmervect)
{
	for(size_t i = 0; i<kmervect.size(); i++)
    {
        kmerdict.insert(make_pair(kmervect[i].rep(), i));
	}
	// cout << "kmervect.rep() in the dict" << endl;
}

int main (int argc, char* argv[]) {

	ifstream filein (argv[1]);
	FILE *fastafile;
	int elem;
	char *buffer;
	std::string kmerstr;
	std::string line;
	Kmer::set_k(KMER_LENGTH);
	dictionary kmerdict;
	std::vector<filedata> allfiles = GetFiles(argv[3]);
    size_t upperlimit = 1000; // 1 thousand reads at a time
    Kmer kmerfromstr;
    Kmers kmervect;
    std::vector<string> seqs;
    std::vector<string> quals;
    Kmers kmersfromreads;
    std::vector<tuple<size_t,size_t,cellspmat>> occurrences;
    std::vector<tuple<size_t,size_t,cellspmat>> transtuples;

	cout << "input file : " << argv[1] <<endl;
	cout << "psbsim depth : 30" << endl;
	cout << "k-mer length : " << KMER_LENGTH <<endl;
	cout << "reference genome : escherichia coli, " << argv[2] <<endl;

	// filtering on kmers --> kmers which occur between 4 and 8 times 
	if(filein.is_open()) {
			while(getline(filein, line)) {
				if(line.length() == 0)
					break;

				string substring = line.substr(1);
				elem = stoi(substring);
				getline(filein, kmerstr); 
				//if(elem>3 && elem<9) {	
				kmerfromstr.set_kmer(kmerstr.c_str());
				kmervect.push_back(kmerfromstr);								
			}
	} else std::cout << "unable to open the input file\n";
	filein.close();

	dictionaryCreation(kmerdict, kmervect);

    size_t read_id = 0; // read_id needs to be global (not just per file)
    for(auto itr=allfiles.begin(); itr!=allfiles.end(); itr++) {

        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);
        
        size_t fillstatus = 1;
        while(fillstatus) {

            fillstatus = pfq->fill_block(seqs, quals, upperlimit);
            size_t nreads = seqs.size();
            
            for(size_t i=0; i<nreads; i++) 
            {
            	// remember that the last valid position is length()-1
                size_t len = seqs[i].length();
                
                // skip this sequence if the length is too short
                if (len <= KMER_LENGTH)
                    continue;

                for(size_t j=0; j<=len-KMER_LENGTH; j++)  
                {

                    std::string kmerstrfromfastq = seqs[i].substr(j, KMER_LENGTH);
                    Kmer mykmer(kmerstrfromfastq.c_str());
                    // remember to use only ::rep() when building kmerdict as well
                    Kmer lexsmall = mykmer.rep();      

                    auto found = kmerdict.find(lexsmall);
                    if(found != kmerdict.end()) {
                    	occurrences.push_back(std::make_tuple(read_id, found->second, make_pair(found->second, vector<size_t>(1,j)))); // vector<tuple<read_id,kmer_id,pair<kmer_id,pos_in_read>>
                        transtuples.push_back(std::make_tuple(found->second, read_id, make_pair(found->second, vector<size_t>(1,j)))); // transtuples.push_back(col_id, row_id, value)
                    }
                }
                read_id++;
            }
        }
    }
    // cout << "fastq file parsed\nsearch ended : vector<tuple<read_id,kmer_id,pair<kmer_id,pos_in_read>> created" << endl;

    CSC<size_t, cellspmat> spmat(occurrences, read_id, kmervect.size(), 
                            [] (cellspmat & c1, cellspmat & c2) 
                            {   if(c1.first != c2.first) cout << "error in MergeDuplicates()" << endl;
                                vector<size_t> merged;
                                merge(c1.second.begin(), c1.second.end(), c2.second.begin(), c2.second.end(), back_inserter(merged));
                                return make_pair(c1.first, merged);
                            });
    std::cout << "spmat created" << endl;

    CSC<size_t, cellspmat> transpmat(transtuples, kmervect.size(), read_id, 
                            [] (cellspmat & c1, cellspmat & c2) 
                            {   if(c1.first != c2.first) cout << "error in MergeDuplicates()" << endl;
                                vector<size_t> merged;
                                merge(c1.second.begin(), c1.second.end(), c2.second.begin(), c2.second.end(), back_inserter(merged));
                                return make_pair(c1.first, merged);
                            });
    std::cout << "transpmat created" << endl;

    // TO DO: define mult and add operation 
   
    spmat.Sorted();
    transpmat.Sorted();
    // std::cout << "spmat and transpmat sorted" << endl;

    std::vector<size_t>	*RowIdsofC = new vector<size_t>[read_id]; // rowidsofc.size() should be correspond to read_id
    std::vector<multcell> *ValuesofC = new vector<multcell>[read_id]; // check valuesofc.size() 

    double start = omp_get_wtime();

    HeapSpGEMM(spmat, transpmat,
        [] (cellspmat & c1, cellspmat & c2) 
        {   if(c1.first != c2.first) cout << "error in multop()" << endl;
            multcell value;
            pair<vector<size_t>, vector<size_t>> vectemp = make_pair(c1.second, c2.second);
            value.insert(make_pair(c1.first, vectemp));
            return value;
        },
        [] (multcell & h, multcell & m)
        {   m.insert(h.begin(), h.end());;
            return m;
        }, 
        RowIdsofC, ValuesofC);
    std::cout << "HeapSpGEMM ok!" << endl;

	return 0;
} 
    /* Ottimizzazione: sfruttare celle diagonali
    Due putnatori che partono dalla cella in alto a sinistra, uno avanza sulla riga e uno avanza sulla colonna 
    Prendo la prima riga */













