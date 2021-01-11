#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <stdio.h>
#include "defs.h"

/* Global Variables */

/* The number of vertices, max. 'clique size', and max.
 * integer 'edge weight' are expressed by default in terms of SCALE.
 * User-defined values can be specified in a config file */
extern LONG_T SCALE;

/* Total number of vertices */
extern LONG_T TotVertices;

/* The generated graph is a collection of highly interconnected
 * 'cliques' with 'inter-clique' edges added probabilistically 
 * This parameter sets the maximum number of vertices in a 'clique' */
extern LONG_T MaxCliqueSize;

/* Probability that two vertices in a 'clique' are connected by
 * a single directed edge, as opposed to edges in both directions.
 * If this value is set to 0, the cluster would become a 'clique'.
 * Self edges are not generated, and so a 'clique' thus defined 
 * with n vertices would have n*(n-1) directed edges */
extern double ProbUnidirectional;

/* initial probability of connecting two 'cliques'
 * inter-clique edges are assigned using a random distribution 
 * that represents a hierarchical thinning as cliques become further apart
 * refer inter-clique edge generation code in graph.c for details */
extern double ProbIntercliqueEdges;

/* maximum number of parallel edges between any two vertices in the graph *
 * note that we generate a multigraph */ 
extern LONG_T MaxParallelEdges;

/* Max. and min. integer edge weights. Edge weights are 
 * randomly chosen from [MinWeight, MaxWeight) */
extern WEIGHT_T MaxWeight;
extern WEIGHT_T MinWeight;

/* Sort edge list by start vertex? */
extern int SORT_EDGELISTS;

/* Sorting alg. to use: 0 for counting sort, 1 for heap sort */
extern int SORT_TYPE;

/* Should the graph be written to file? */	
extern int WRITE_TO_FILE;

/* Default output file */ 
extern char OUTFILE[30];

/* Default log file, for printing auxiliary information */
extern char LOGFILE[30];

#endif
