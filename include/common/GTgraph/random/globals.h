#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <stdio.h>
#include "defs.h"

/* Global Variables */

/* The random graph model to use
 * GRAPH_MODEL 0 : G(n, m) -- m edges are added choosing
 *                 n pairs of vertices at random  
 * GRAPH_MODEL 1 : G(n, p) -- Erdos-Renyi graph  */
extern int GRAPH_MODEL;

/* Total number of vertices */
extern LONG_T n;

/* m, the no. of edges in graph model 0 */
extern LONG_T m;

/* Probability of an edge between any two vertices in an 
 * Erdos-Renyi random graph */
extern double p;

/* Max. and min. integer edge weights. Edge weights are 
 * randomly chosen from [MIN_WEIGHT, MAX_WEIGHT) */
extern WEIGHT_T MAX_WEIGHT;
extern WEIGHT_T MIN_WEIGHT;

/* Are self loops allowed? */
extern int SELF_LOOPS;

/* Sort edge list by start vertex before writing to file */
extern int SORT_EDGELISTS;

/* Sorting alg. to use: 0 for counting sort, 1 for heap sort */
extern int SORT_TYPE;

/* Should the graph be written to file? */	
extern int WRITE_TO_FILE;

/* If this option is selected, memory is allocated to the 
   graph data structure in addition to / instead of 
   writing the graph to disk */
extern int STORE_IN_MEMORY;

/* Default output file */ 
extern char OUTFILE[30];

/* Default log file, for printing auxiliary information */
extern char LOGFILE[30];

#endif
