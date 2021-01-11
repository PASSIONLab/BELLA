#ifndef _GRAPH_H
#define _GRAPH_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "sprng.h"
#include "globals.h"
#include "defs.h"
#include "utils.h"

#define SPRNG_SEED1 12619830
#define SPRNG_SEED2 31219885
#define SPRNG_SEED3 72824922
#define SPRNG_SEED4 81984016

#define DEBUG 0

extern int *stream1, *stream2, *stream3, *stream4;

void graphGen(graph*);
void choosePartition(LONG_T*, LONG_T*, LONG_T);
void varyParams(double*, double*, double*, double*);
#endif
