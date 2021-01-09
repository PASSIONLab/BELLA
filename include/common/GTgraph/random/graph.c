#include "graph.h"

void graphGen(graph* g)
{

	LONG_T *startVertex, *endVertex;
	WEIGHT_T *weights;
	LONG_T i, j, u, v;
	WEIGHT_T w;
	LONG_T estNumEdges, numEdges, edgeNum;	
	int *stream1, *stream2;
	FILE* outfp;

	/*----------------------------------------------*/
	/*		initialize SPRNG 		*/
	/*----------------------------------------------*/

	stream1 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
	stream2 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED2, SPRNG_DEFAULT);

	/*------------------------------------------------------*/
	/*		generate edges as per the		*/		 
	/*		graph model and user options	      	*/
	/*------------------------------------------------------*/

	if ((STORE_IN_MEMORY == 0) && (SORT_EDGELISTS == 0)) {
		fprintf(stderr, "Generating edges on the fly\n");
		outfp = fopen(OUTFILE, "w");
                fprintf(outfp, "c FILE			: %s\n", OUTFILE);
                fprintf(outfp, "c No. of vertices	: %ld\n", n);
		if (GRAPH_MODEL == 1) 
			fprintf(outfp, "c No. of edges		: %ld\n", m);
		else
			fprintf(outfp, "                                                 \n");
                fprintf(outfp, "c Max. weight	        : %ld\n", MAX_WEIGHT);
                fprintf(outfp, "c Min. weight          	: %ld\n", MIN_WEIGHT);
                fprintf(outfp, "c A directed arc from u to v of weight w\n");
                fprintf(outfp, "c is represented below as ' a  u  v  w '\n");
		fprintf(stderr, "Generating Edges ... ");
		
		/* Erdos-Renyi */ 	
		if (GRAPH_MODEL == 0) {
			
			/* Write the no. of edges later */
			fprintf(outfp, "                             \n");
			numEdges = 0;
			
			for (i=0; i<n; i++) {
				for (j=0; j<n; j++) {
					
					if ((i==j) && (SELF_LOOPS == 0))		
						continue;

					if (p > sprng(stream1)) {
						w = MIN_WEIGHT + (WEIGHT_T) \
						   (MAX_WEIGHT - MIN_WEIGHT) * sprng(stream2);
						/* Create edge */
						fprintf(outfp, "a %ld %ld %ld\n", i+1, j+1, w);
						numEdges++;
					}
							
				} 
			}

			m = numEdges;	
			fclose(outfp);
			fprintf(stderr, "done\n");
			updateEdgeVal(OUTFILE, numEdges);			

		} else {

                        fprintf(outfp, "p sp %ld %ld\n", n, m);
                        numEdges = 0;

                        for (i=0; i<m; i++) {

				u = (LONG_T) isprng(stream1) % n;
				v = (LONG_T) isprng(stream1) % n;
				if ((u == v) && (SELF_LOOPS == 0)) {
					i--;
					continue;
				}

				w = MIN_WEIGHT + (WEIGHT_T) (MAX_WEIGHT - MIN_WEIGHT) * sprng(stream2);

                                /* Create edge */
                                fprintf(outfp, "a %ld %ld %ld\n", u+1, v+1, w);
                               
                        }

			fclose(outfp);
                        fprintf(stderr, "done\n");
	
		}	

		free(stream1);
		free(stream2);
		return;

	}

	fprintf(stderr, "Generating edges ... ");

	if (GRAPH_MODEL == 0) {

		/* Estimate the no. of edges */
		if (SELF_LOOPS)
			estNumEdges = (LONG_T) (120 * n * n * p)/100;
		else
			estNumEdges = (LONG_T) (120 * n * (n-1) *  p)/100; 

		edgeNum = 0;
		numEdges = 0;

		startVertex = (LONG_T *) malloc(estNumEdges * sizeof(LONG_T)); 
		endVertex   = (LONG_T *) malloc(estNumEdges * sizeof(LONG_T));
		
		for (i=0; i<n; i++) {
			
			for (j=0; j<n; j++) {

				if ((i==j) && (SELF_LOOPS == 0))
                                	continue;

                                if (p > sprng(stream1)) {

					startVertex[edgeNum] = i; 
					endVertex[edgeNum]   = j;
					edgeNum++;

				}	
			}

		}

		numEdges = edgeNum;

	} else {

		startVertex = (LONG_T *) malloc(m * sizeof(LONG_T));
		endVertex   = (LONG_T *) malloc(m * sizeof(LONG_T));

		for (i=0; i<m; i++) {

			u = (LONG_T) isprng(stream1) % n;
               		v = (LONG_T) isprng(stream1) % n;
                	if ((u == v) && (SELF_LOOPS == 0)) {
				i--;
                		continue;
			}
			startVertex[i] = u;
			endVertex[i] = v;
		}		
		
		numEdges = m;	

	}

	fprintf(stderr, "done\n");

	free(stream1);

	/*----------------------------------------------*/
        /*              generate edge weights           */
        /*----------------------------------------------*/

        fprintf(stderr, "Generating edge weights ... ");

        weights = (WEIGHT_T *) malloc(numEdges*sizeof(WEIGHT_T));

        for (i=0; i<numEdges; i++) {
                weights[i] = MIN_WEIGHT + (WEIGHT_T) (MAX_WEIGHT - MIN_WEIGHT) \
                        * sprng(stream2);
        }

        fprintf(stderr, "done\n");
        free(stream2);

        /*-------------------------------------------------------*/
        /*              sort the edge lists with start           */
        /*              vertex as primary key                    */
        /*---------------------------------------------- --------*/

        if (SORT_EDGELISTS) {

                fprintf(stderr, "Sorting edge list by start vertex ... ");
		if (GRAPH_MODEL == 1) {
                	if (SORT_TYPE == 0) {
                        	/* Counting sort */
                        	countingSort(startVertex, endVertex, weights, numEdges);
                	} else {
                        	/* Heap sort */
                       		heapSort(startVertex, endVertex, weights, numEdges);
                	}
		}
                fprintf(stderr, "done\n");
        }

        g->start = startVertex;
        g->end = endVertex;
        g->w = weights;
        g->n = n;
        g->m = numEdges;
}
