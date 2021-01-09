#include <stdio.h>
#include <stdlib.h>
#include "graph.h"
#include "init.h"

int main(int argc, char** argv)
{
	
	/* data structure for representing the generated edge sets 
	 * see defs.h */
	graph G;

	/* init.c */
 	parseUserInput(argc, argv);
	
	/* graph.c */
	graphGen(&G);

	/* utils.c */
	writeToFile(&G);

	/* Free memory */
	if (STORE_IN_MEMORY) {
		free(G.start);
		free(G.end);
		free(G.w);
	}

	return(0);	
}
