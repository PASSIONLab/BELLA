#ifndef _INIT_H
#define _INIT_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "defs.h"
#include "utils.h"
#include "globals.h"
#include "string.h"

void parseUserInput(int, char**);
void usage(void);

void getParamsFromFile(char*);
void getParams(void);
void setGTgraphParams(int scale, int edgefactor, double a0, double b0, double c0, double d0);

#endif

