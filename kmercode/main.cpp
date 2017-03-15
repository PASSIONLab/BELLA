#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>
using namespace std;

#define KMERLONGS (MAX_KMER_SIZE/32)      // 32 = numbits(uint64_t)/2-  with 2 being the number of bits needed per nucleotide
#define KMER_LENGTH 15

#include "Kmer.hpp"
#include "common.h"
#include "ParallelFASTQ.h"

struct filedata
{
    char filename[MAX_FILE_PATH];
    size_t filesize;
};



vector< filedata >  GetFiles(char * filename)
{
    int64_t totalsize = 0;
    int numfiles = 0;
    vector< filedata > filesview;
    
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
        cout << filesview.back().filename << " : " << filesview.back().filesize / (1024*1024) << " MB" << endl;
        allfiles.getline(fdata.filename,MAX_FILE_PATH);
        totalsize += fdata.filesize;
        numfiles++;
    }
    return filesview;
}

int main(int argc, char* argv[])
{
    Kmer::set_k(KMER_LENGTH);
    vector<filedata> allfiles = GetFiles(argv[1]);
    size_t upperlimit = 1000;   // a thousand reads at a time
    
    vector<string> seqs;
    vector<string> quals;
    for(auto itr= allfiles.begin(); itr != allfiles.end(); itr++)
    {
        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);
        
        size_t fillstatus = 1;
        while(fillstatus)
        {
            fillstatus = pfq->fill_block(seqs, quals, upperlimit);
            size_t nreads = seqs.size();
            
            for(size_t i=0; i< nreads; ++i)
            {
                size_t found;
                found = seqs[i].length();
                
                // skip this sequence if the length is too short
                if (seqs[i].length() <= KMER_LENGTH) {
                    continue;
                }
                int nkmers = (seqs[i].length()-KMER_LENGTH+1);

                
                std::vector<Kmer> kmers = Kmer::getKmers(seqs[i]); // calculate all the kmers
                assert(kmers.size() == nkmers);
                
                // Aydin: the following lines are just for debugging, comment out after you get a feeling of how it works
                cout << "printing the first k-mer and its reverse complement from this read: ";
                cout << kmers[0].toString() << " " << kmers[0].twin().toString() << endl;
            }
            
        }
    }
}
