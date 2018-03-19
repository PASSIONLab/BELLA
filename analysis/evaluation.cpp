#include "evaluation.h" 
#include "IntervalTree.h"
#include "../mtspgemm2017/global.h"
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

    //
    // Program name and purpose
    //
    cout << "\nProgram to test BELLA (and others) sensitivy and precision" << endl;

    //
    // Setup the input files
    //
    option_t *optList, *thisOpt;
    // Get list of command line options and their arguments 
    optList = NULL;
    optList = GetOptList(argc, argv, (char*)"g:b:m:p:d:hl:z");
   
    int ovl_len = 2000;
    bool simulated = false;
    char *truth = NULL;
    char *bellafile = NULL;
    char *minimapfile = NULL;   
    char *mhapfile = NULL;
    char *dalignerfile = NULL; 
    char *blasrfile = NULL; 
  
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
                truth = strdup(thisOpt->argument);
                break;
            }
            case 'b': {
                if(thisOpt->argument == NULL)
                {
                    cout << "\nProgram execution terminated: -b requires BELLA output file" << endl;
                    cout << "Run with -h to print out the command line options\n" << endl;
                    return 0;
                }
                bellafile = strdup(thisOpt->argument);
                break;
            }
            case 'z': simulated = true; break; // Using simulated reads
            case 't': {
                ovl_len = atoi(thisOpt->argument);
                break;
            }
            case 'm': {
                minimapfile = strdup(thisOpt->argument);
                break;
            }
            case 'p': {
                mhapfile = strdup(thisOpt->argument);
                break;
            }
            case 'd': {
                dalignerfile = strdup(thisOpt->argument);
                break;
            }
            case 'l': {
                blasrfile = strdup(thisOpt->argument);
                break;
            }
            case 'h': {
                cout << "\nUsage:\n" << endl;
                cout << " -z : Simulated reads [false]" << endl;
                cout << " -g : ground truth file (required)" << endl;
                cout << " -n : number of reads of fastq(s) (required)" << endl; // need to have this value to 
                                                                                // compute minimap sensitivity (impractical but 
                                                                                // safer as with real data isn't safe counting 
                                                                                // the number of reads from the ground truth)
                cout << " -t : Overlap length [2000]" << endl;
                cout << " -b : BELLA output file (required)" << endl;
                cout << " -m : Minimap or Minimap2 output file" << endl;
                cout << " -p : MHAP output file" << endl;
                cout << " -d : Daligner output file" << endl;
                cout << " -l : Blasr output file\n" << endl;
                
                FreeOptList(thisOpt); // Done with this list, free it
                return 0;
            }
        }
    }

    if(truth == NULL || bellafile == NULL)
    {
        cout << "\nProgram execution terminated: missing arguments" << endl;
        cout << "Run with -h to print out the command line options\n" << endl;
        return 0;
    }

    free(optList);
    free(thisOpt);

    std::ifstream ground(truth);
	std::ifstream bella(bellafile);
    std::ifstream minimap(minimapfile);
    std::ifstream mhap(mhapfile);
    std::ifstream blasr(blasrfile);
    std::ifstream daligner(dalignerfile);

    benchmarkingAl(ground, bella, minimap, mhap, blasr, daligner, simulated, ovl_len);  

}