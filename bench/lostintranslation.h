#include "utils.h"
#include "common.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <sys/types.h> 
#include <sys/stat.h> 
#include <math.h>
#include <limits.h>
#include <bitset>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h> 
#include <sstream>
#include <set>
#include <memory>
#include <typeinfo>
#include <omp.h>

using namespace std;

/* PAF format: https://github.com/lh3/miniasm/blob/master/PAF.md
 * column seq name
 * column seq length
 * column seq start
 * column seq end
 * strand (+/-)
 * row seq name
 * row seq length
 * row seq start
 * row seq end
 * number of residue matches (alignment score)  ---> GGGG: I set this to 0 if missing
 * alignment block length (overlap length)      ---> GGGG: I compute this if missing
 * mapping quality (0-255; 255 for missing) */

/* If PAF is generated from an alignment, column 10 equals the number of sequence matches, and column 11 equals the total 
number of sequence matches, mismatches, insertions and deletions in the alignment. If alignment is not available, 
column 10 and 11 are still required but may be highly inaccurate. */

int estimate (int begpV, int endpV, int lenV, int begpH, int endpH, int lenH)
{
    int diffV = begpV - begpV;
    int diffH = endpH - begpH;
    int minL  = min(begpV, begpH);
    int minR  = min(lenV - endpV, lenH - endpH);
    int ovlen = minL + minR + (diffV + diffH)/2;

    return ovlen;
}

vector<std::string> split (const std::string &s, char delim)
{
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (std::getline (ss, item, delim))
    {
        result.push_back (item);
    }

    return result;
}

void MHAP2PAF(ifstream& input, char* filename)
{
    int maxt = 1;
#pragma omp parallel
    {
        maxt = omp_get_num_threads();
    }

    uint64_t numoverlap = std::count(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>(), '\n');
    input.seekg(0, std::ios_base::beg);

    vector<std::string> entries;
    vector<std::stringstream> local(maxt);      

    /* read input file */
    if(input)
        for (int i = 0; i < numoverlap; ++i)
        {
            std::string line;
            std::getline(input, line);
            entries.push_back(line);
        }
    input.close();

    /* transform MHAP output in PAF format */
#pragma omp parallel for
    for(uint64_t i = 0; i < numoverlap; i++) 
    {
        std::stringstream linestream(entries[i]);
        int ithread = omp_get_thread_num();

        /* MHAP format: cname, rname, err, nkmer, cstrand, cstart, cend, clen, rstrand, rstart, rend, rlen */
        std::vector<std::string> v = split (entries[i], '\t');

        /* change strand formatting */
        if(v[8] == "0") v[8] = "+";         
            else v[8] = "-";

        /* compute overlap length if missing (begpV, endpV, lenV, begpH, endpH, lenH) */
        int ovlen = estimate (stoi(v[5]), stoi(v[6]), stoi(v[7]), stoi(v[9]), stoi(v[10]), stoi(v[11]));

        /* GGGG: improve readability */
        local[ithread] << v[0] << "\t" << v[7] << "\t" << v[5] << "\t" << v[6] << "\t" << v[8] 
            << "\t" << v[1] << "\t" << v[11] << "\t" << v[9] << "\t" << v[10] << "\t" << "0" 
                << "\t" << ovlen << "\t255" << endl;
    }

    /* write to a new file */
    int64_t * bytes = new int64_t[maxt];
    for(int i = 0; i < maxt; ++i)
    {
        local[i].seekg(0, ios::end);
        bytes[i] = local[i].tellg();
        local[i].seekg(0, ios::beg);
    }
    int64_t bytestotal = std::accumulate(bytes, bytes + maxt, static_cast<int64_t>(0));

    std::ofstream output(filename, std::ios::binary | std::ios::app);
#ifdef PRINT
    cout << "Creating or appending to output file with " << (double)bytestotal/(double)(1024 * 1024) << " MB" << endl;
#endif
    output.seekp(bytestotal - 1);
    /* this will likely create a sparse file so the actual disks won't spin yet */
    output.write("", 1); 
    output.close();

    #pragma omp parallel
    {
        int ithread = omp_get_thread_num(); 

        FILE *ffinal;
        /* then everyone fills it */
        if ((ffinal = fopen(filename, "rb+")) == NULL) 
        {
            fprintf(stderr, "File %s failed to open at thread %d\n", filename, ithread);
        }
        int64_t bytesuntil = std::accumulate(bytes, bytes + ithread, static_cast<int64_t>(0));
        fseek (ffinal , bytesuntil , SEEK_SET);
        std::string text = local[ithread].str();
        fwrite(text.c_str(),1, bytes[ithread], ffinal);
        fflush(ffinal);
        fclose(ffinal);
    }
    delete [] bytes;
}