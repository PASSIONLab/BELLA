/**
* Evaluation code to compare long read-to-read aligner recall/precision
* @author: Giulia Guidi
*/
#ifdef __cplusplus
extern "C" {
#endif

#include "../optlist/optlist.h" // command line parser

#ifdef __cplusplus
}
#endif

#include "truth.h"
#include "overlapper.h" 
#include "onlyoverlap.h" 
#include "IntervalTree.h"
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
#include <string.h>
#include <unistd.h>
#include <ctype.h> 
#include <memory>

int main (int argc, char* argv[]) {

    // Program name and purpose
    cout << "\nProgram to Compute Long-Read Overlapper Recall, Precision, and Specificity" << endl;
    // Setup the input files
    option_t *optList, *thisOpt;
    // Get list of command line options and their arguments 
    optList = NULL;
    optList = GetOptList(argc, argv, (char*)"g:b:a:m:p:d:hl:zo:c:i:O");

    int ovLen = 2000;   // min overlap length to be considered a true positive
    bool sim = false;   // simulated dataset [false]
    bool oov = false;
    char *th = NULL;    // truth
    char *b = NULL;     // bella
    char *a = NULL;     // BELLA in PAF format
    char *o = NULL;     // BELLA evaluation output filename
    char *m = NULL;     // minimap/miniamap2
    char *p = NULL;     // mhap
    char *d = NULL;     // daligner
    char *l = NULL;     // blasr
    char *c = NULL;     // mecat
    char *i = NULL;     // MECAT idx2read filename

    if(optList == NULL)
    {
        cout << "Program execution terminated: not enough parameters or invalid option" << endl;
        cout << "Run with -h to print out the command line options" << endl;
        return 0;
    }

    while (optList!=NULL) {
        thisOpt = optList;
        optList = optList->next;
        switch (thisOpt->option) {
            case 'g': {
                if(thisOpt->argument == NULL)
                {
                    cout << "\nProgram execution terminated: -g requires the ground truth file" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                th = strdup(thisOpt->argument);
                break;
            }
            case 'b': {
                b = strdup(thisOpt->argument);
                break;
            }
            case 'a': {
                a = strdup(thisOpt->argument);
                break;
            }
            case 'z': sim = true; break; // using simulated data
            case 'O': oov = true; break; // using simulated data
            case 't': {
                ovLen = atoi(thisOpt->argument);
                break;
            }
            case 'm': {
                m = strdup(thisOpt->argument);
                break;
            }
            case 'p': {
                p = strdup(thisOpt->argument);
                break;
            }
            case 'd': {
                d = strdup(thisOpt->argument);
                break;
            }
            case 'l': {
                l = strdup(thisOpt->argument);
                break;
            }
            case 'o': {
                o = strdup(thisOpt->argument);
                break;
            }
            case 'c': {
                c = strdup(thisOpt->argument);
                break;
            }
            case 'i': {
                i = strdup(thisOpt->argument);
                break;
            }
            case 'h': {
                cout << "\nUsage:\n" << endl;
                cout << " -z : Simulated data [false]" << endl;
                cout << " -g : Ground truth file (required)" << endl;
                cout << " -t : Overlap length [2000]" << endl;
                cout << " -O : BELLA output only overlap [false]" << endl;
                cout << " -b : BELLA output file" << endl;
                cout << " -a : BELLA output file in PAF format" << endl;
                cout << " -m : MINIMAP2 output file" << endl;
                cout << " -p : MHAP output file" << endl;
                cout << " -d : DALIGNER output file" << endl;
                cout << " -l : BLASR output file\n" << endl;
                cout << " -t : MECAT output file\n" << endl;
                cout << " -i : MECAT idx2read file\n" << endl;
                cout << " -o : Filename where to output extra information from BELLA's evaluation\n" << endl;
                FreeOptList(thisOpt); // Done with this list, free it
                return 0;
            }
        }
    }

    if(th == NULL)
    {
        cout << "\nProgram execution terminated: missing argument." << endl;
        cout << "Run with -h to print out the command line options.\n" << endl;
        return 0;
    }

    if(b != NULL && o == NULL)
    {
        cout << "\nProgram execution terminated: missing argument." << endl;
        cout << "Please add name for output file with -o option. Run with -h to print out the command line options.\n" << endl;
        return 0;
    }

    if(a != NULL && o == NULL)
    {
        cout << "\nProgram execution terminated: missing argument." << endl;
        cout << "Please add name for output file with -o option. Run with -h to print out the command line options.\n" << endl;
        return 0;
    }

    if(c != NULL && i == NULL)
    {
        cout << "\nProgram execution terminated: missing argument." << endl;
        cout << "Please add MECAT idx2read file with -i option. Run with -h to print out the command line options.\n" << endl;
        return 0;
    }

    cout << "\nDefinitions:" << endl;
    cout << "   .Recall = TP/(TP+FN)"  << endl;
    cout << "   .Precision = TP/(TP+FP)" << endl;
    cout << "   .Specificity = TN/(TN+FP)" << endl;

    free(optList);
    free(thisOpt);

    std::ifstream thf(th);
    multimap<string,readInfo> seqmap;
    uint32_t num = getTruth(thf,sim,ovLen,seqmap,num); // will contain ground truth cardinality

    if(oov)
    {
        if(b != NULL)
        {
            std::ifstream bf(b);
            myBella(thf,bf,sim,ovLen,seqmap,num); // bella only overlaps (not alignment score, etc.)
        }
    }
    else
    {
        if(b != NULL)
        {
            std::ifstream bf(b);
            std::string out(o);
            metricsBella(thf,bf,sim,ovLen,out,seqmap,num); // bella
        }
    }
    if(a != NULL)
    {
        std::ifstream bf(a);
        std::string out(o);
        metricsBellaPAF(thf,bf,sim,ovLen,out,seqmap,num); // bella in paf format
    }
    if(m != NULL)
    {
        std::ifstream mf(m);
        metricsMinimap2(thf,mf,sim,ovLen,seqmap,num); // minimap/minimap2
    }
    if(p != NULL)
    {
        std::ifstream pf(p);
        metricsMhap(thf,pf,sim,ovLen,seqmap,num); // mhap
    }
    if(l != NULL)
    {
        std::ifstream lf(l);
        metricsBlasr(thf,lf,sim,ovLen,seqmap,num); // blasr
    }
    if(d != NULL)
    {
        std::ifstream df(d);
        metricsDaligner(thf,df,sim,ovLen,seqmap,num); // daligner
    }
    if(c != NULL)
    {
        std::ifstream cf(c);
        std::ifstream idx2read(i);
        metricsMecat(thf,cf,sim,ovLen,seqmap,num,idx2read); // mecat
    }
}
