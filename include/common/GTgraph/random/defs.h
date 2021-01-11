#ifndef _DEFS_H
#define _DEFS_H

/* Using longs to represent vertex numbers */
#define LONG_T long int

/* Edge Weights are also longs */
#define WEIGHT_T long int

typedef struct 
{
	/* No. of edges, represented by m */
	LONG_T m;

	/* No. of vertices, represented by n */
	LONG_T n;
	
	/* Arrays of size 'm' storing the edge information
	 * A directed edge 'e' (0 <= e < m) from start[e] to end[e]
	 * had an integer weight w[e] */
	LONG_T* start;
	LONG_T* end;
	WEIGHT_T* w;

} graph;

#endif
