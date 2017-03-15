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
    
    for(auto itr= allfiles.begin(); itr != allfiles.end(); itr++)
    {
        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(itr->filename, false, itr->filesize);
    }
}
