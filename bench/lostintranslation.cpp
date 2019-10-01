
//=======================================================================
// Title:  C++ program to translate overlappers' output in PAF format
// Author: G. Guidi
// Date:   23 Feb 2019
//=======================================================================

#ifdef __cplusplus
extern "C" {
#endif
#include "../optlist/optlist.h" /* command line parser */
#ifdef __cplusplus
}
#endif

#include "lostintranslation.h"  
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

using namespace std;

int main (int argc, char* argv[]) {

    /* program description */
    cout << "\nProgram to translate overlappers' outputs in PAF format\nRun with a single overlapper's output (for now)" << endl;
    /* input files setup */
    option_t *optList, *thisOpt;
    /* list of command line options and their arguments */
    optList = NULL;
    optList = GetOptList(argc, argv, (char*)"b:B:m:t:r:d:f:i:h");

    char *b = NULL;         // bella
    // char *B = NULL;      // dibella // not supported yet
    char *m = NULL;         // mhap
    char *t = NULL;         // mecat
    char *r = NULL;         // blasr
    char *d = NULL;         // daligner
    char *filename = NULL;  // filename translated output
    char *index = NULL;     // filename translated output

    if(optList == NULL)
    {
        cout << "Program execution terminated: not enough parameters or invalid option" << endl;
        cout << "Run with -h to print out the command line options" << endl;
        return 0;
    }

    while (optList!=NULL)
    {
        thisOpt = optList;
        optList = optList->next;
        switch (thisOpt->option)
        {
            case 'b': {
                b = strdup(thisOpt->argument);
                break;
            }
            //case 'B': {
            //    B = strdup(thisOpt->argument);
            //    break;
            //}
            case 'm': {
                m = strdup(thisOpt->argument);
                break;
            }
            case 't': {
                t = strdup(thisOpt->argument);
                break;
            }
            case 'i': {
                index = strdup(thisOpt->argument);
                break;
            }
            case 'r': {
                r = strdup(thisOpt->argument);
                break;
            }
            case 'd': {
                d = strdup(thisOpt->argument);
                break;
            }
            case 'f': {
                filename = strdup(thisOpt->argument);
                break;
            }
            case 'h': {
                cout << "\nUsage:\n" << endl;
                cout << " -b : BELLA" << endl;
                // cout << " -B : diBELLA" << endl;
                cout << " -m : MHAP" << endl;
                cout << " -t : MECAT" << endl;
                cout << " -i : MECAT index" << endl;
                cout << " -r : BLASR" << endl;
                cout << " -d : DALIGNER" << endl;
                cout << " -f : filename" << endl;
                cout << " -h : usage\n" << endl;
                /* done with this list, free it */
                FreeOptList(thisOpt);
                return 0;
            }
        }
    }

    free(optList);
    free(thisOpt);

    if(b != NULL)
    {
        ifstream input(b);
        BELLA2PAF(input, filename);
    }
    else if(m != NULL)
    {
        ifstream input(m);
        MHAP2PAF(input, filename);
    }
    else if(r != NULL)
    {
        ifstream input(r);
        BLASR2PAF(input, filename);
    }
    else if(t != NULL)
    {
        if(index != NULL)
        {
            ifstream input(t);
            ifstream idx(index);
            MECAT2PAF(input, filename, idx);
        }
        else
        {
            std::cout << "MECAT index file is missing" << endl;
            return 0;
        }
    }
    else if(d != NULL)
    {
        ifstream input(d);
        DALIGNER2PAF(input, filename);
    }
    //else if(B != NULL)
    //{
    //    ifstream input(B);
    //    diBELLA2PAF(input, filename);
    //}
    return 0;
}
