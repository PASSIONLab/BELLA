#ifndef _UTILS_H
#define _UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "defs.h"
#include "globals.h"
#include "string.h"

#define bits(x,k,j) ((x>>k) & ~(~0<<j))

void writeToFile(graph*);
void updateLog(void);
void countingSort(LONG_T*, LONG_T*, WEIGHT_T*, LONG_T);
void countsort_aux(LONG_T q, LONG_T *lKey, LONG_T *lSorted, 
		   LONG_T * auxKey, LONG_T* auxSorted,
		   LONG_T R, LONG_T bitOff, LONG_T m);

/* heap sort code adapted from 
 * http://linux.wku.edu/~lamonml/algor/sort/heap.html */
void heapSort(LONG_T*, LONG_T*, WEIGHT_T*,  LONG_T);
void siftDown(LONG_T*, LONG_T*, WEIGHT_T*, LONG_T, LONG_T);
void swapL(LONG_T*, LONG_T*);
void swapW(WEIGHT_T*, WEIGHT_T*);
#endif
