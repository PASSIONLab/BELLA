#include "graph.h"

void graphGen(graph* g)
{

	LONG_T* clusterSizes;
	LONG_T* firstVsInCluster;
	LONG_T estTotClusters, totClusters;
	
	LONG_T *startVertex, *endVertex;
	LONG_T numEdges;
        LONG_T numIntraClusterEdges, numInterClusterEdges;
	WEIGHT_T* weights;
	LONG_T i_cluster, currCluster;
	LONG_T *startV, *endV, *d;
	LONG_T estNumEdges, edgeNum;

	LONG_T i, j, k, t, t1, t2, dsize;
	double p;
	FILE* outfp;
	int* stream;
	LONG_T* permV;


	/*----------------------------------------------*/
	/*		initialize SPRNG 		*/
	/*----------------------------------------------*/

	stream = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED, SPRNG_DEFAULT);

	/*----------------------------------------------*/
	/*		generate clusters		*/
	/*----------------------------------------------*/

  	/* Estimate number of clusters required to make up 
	 * TotVertices and pad by 25% */
	estTotClusters = 1.25 * TotVertices / (MaxCliqueSize/2);
	
	/* Allocate a block of memory for this cluster-size array */
	clusterSizes = (LONG_T *) malloc(estTotClusters*sizeof(LONG_T));

	/* Generate random cluster sizes. */
	for(i = 0; i < estTotClusters; i++) {
		clusterSizes[i] = 1 + (sprng(stream) * MaxCliqueSize);
	}

	totClusters = 0;
  
	/* Allocate memory for storing the first vertex in each cluster */
	firstVsInCluster = (LONG_T *) malloc(estTotClusters*sizeof(LONG_T));

	/* Compute the first vertex in each cluster */
	firstVsInCluster[0] = 0;
	for (i=1; i<estTotClusters; i++) {
		firstVsInCluster[i] = firstVsInCluster[i-1] + clusterSizes[i-1];
		if (firstVsInCluster[i] > TotVertices-1)
			break;
	}

	totClusters = i;

	/* Fix the size of the last cluster */
	clusterSizes[totClusters-1] = TotVertices -\
		firstVsInCluster[totClusters-1];

	/*------------------------------------------------------*/
	/*		generate intra-cluster edges		*/
	/*------------------------------------------------------*/
 
	/* Roughly estimate the total number of edges */
	estNumEdges = (LONG_T) ((TotVertices * (double) MaxCliqueSize * (2-ProbUnidirectional)/2) +\
		      	        (TotVertices * (double) ProbIntercliqueEdges/(1-ProbIntercliqueEdges))) * \
				(1+MaxParallelEdges/2);

	fprintf(stderr, "\n\tEstimated number of edges - %ld\n", estNumEdges);
	fprintf(stderr, "Generating intra-cluster edges ... ");

	/* Check if no. of edges is within bounds for 32-bit code */
	if ((estNumEdges > ((1<<30) - 1)) && (sizeof(LONG_T*) < 8)) {
		fprintf(stderr, "ERROR: long* should be 8 bytes \
				for this problem size\n");
		fprintf(stderr, "\tPlease recompile the code \
			       	in 64-bit mode\n");
		exit(-1);
	}

	edgeNum = 0;
	p = ProbUnidirectional;
#if DEBUG
	fprintf (stderr, "[allocating %3.3f GB memory ... ", (double) 2*estNumEdges*8/(1<<30));
#endif
	startV = (LONG_T *) malloc(estNumEdges*sizeof(LONG_T));
	endV = (LONG_T *) malloc(estNumEdges*sizeof(LONG_T));
#if DEBUG
	fprintf(stderr, "done] ");  
#endif
	for (i_cluster=0; i_cluster < totClusters; i_cluster++) {

		for (i = 0; i < clusterSizes[i_cluster]; i++) {

			for (j = 0; j < i; j++) {
			
				if (sprng(stream) >= p) {
					/* generate edges in both directions */
					for (k=0; k<1+((LONG_T) (MaxParallelEdges-1) * sprng(stream)); k++) {
						startV[edgeNum] = i + \
							firstVsInCluster[i_cluster];
						endV[edgeNum] = j + \
							firstVsInCluster[i_cluster];
						edgeNum++;
					}
					for (k=0; k<1+((LONG_T) (MaxParallelEdges-1) * sprng(stream)); k++) {
						startV[edgeNum] = j + \
							firstVsInCluster[i_cluster];
                                        	endV[edgeNum] = i + \
							firstVsInCluster[i_cluster];
						edgeNum++;
					}
				} else {

					/* generate an edge only in 
					 * one direction */
					if (sprng(stream) > 0.5) {
						for (k=0; k<1+((LONG_T) (MaxParallelEdges-1) * sprng(stream)); k++) {	
							startV[edgeNum] = i + \
							firstVsInCluster[i_cluster];
							endV[edgeNum] = j + \
							firstVsInCluster[i_cluster];
							edgeNum++;
						}
					} else {
						for (k=0; k<1+((LONG_T) (MaxParallelEdges-1) * sprng(stream)); k++) {
							startV[edgeNum] = j + \
							firstVsInCluster[i_cluster];	
							endV[edgeNum] = i + \
							firstVsInCluster[i_cluster];
							edgeNum++;
						}
					}	
				}
			}
			
		}
	}

	numIntraClusterEdges = edgeNum;

	fprintf(stderr, "done\n");
	printf("\tNo. of intra-cluster edges - %ld \n", numIntraClusterEdges);
	
	/*----------------------------------------------*/
	/*		connect the clusters		*/
	/*----------------------------------------------*/

        fprintf(stderr, "Generating inter-clique edges ... ");	

	/* generate exponential distances as suggested in the spec */
	dsize = (LONG_T) (log((double)TotVertices)/log(2));
	d = (LONG_T *) malloc(dsize * sizeof(LONG_T));
	for (i = 0; i < dsize; i++) {
		d[i] = (LONG_T) pow(2, (double) i);
	}

	currCluster = 0;

	for (i = 0; i < TotVertices; i++) {

		p = ProbIntercliqueEdges;	

		/* determine current cluster */
		for (j = currCluster; j<totClusters; j++) {
			if ((i >= firstVsInCluster[j]) && \
			    (i < firstVsInCluster[j] + clusterSizes[j])) {
				currCluster = j;
				break;
			}	
		}

		for (t = 1; t < dsize; t++) {

			j = (i + d[t] + (LONG_T) (sprng(stream) * (d[t] - d[t-1]))) % TotVertices;	

			/* Ensure that i and j don't belong to the same cluster */
	
			if ((j<firstVsInCluster[currCluster]) || \
			    (j>=firstVsInCluster[currCluster] + clusterSizes[currCluster])) {
				
				for (k=0; k<1+((LONG_T) (MaxParallelEdges - 1)* sprng(stream)); k++) {
					if (p > sprng(stream)) {
						startV[edgeNum] = i;
						endV[edgeNum] = j;
						edgeNum++;	
					}	
				}	
			}
			
			p = p/2;
		}
	}
	
	numEdges = edgeNum;
	numInterClusterEdges = numEdges - numIntraClusterEdges;	

	free(clusterSizes);  
	free(firstVsInCluster);
	free(d);

	fprintf(stderr, "done\n");
	fprintf(stderr, "\tNo. of inter-cluster edges - %ld\n", numInterClusterEdges);
	fprintf(stderr, "\tTotal no. of edges - %ld\n", numEdges);

	/*--------------------------------------------------------------*/
	/*		shuffle vertices to remove locality		*/
	/*--------------------------------------------------------------*/

	fprintf(stderr, "Shuffling vertices to remove locality ... ");
#if DEBUG
	fprintf (stderr, "[allocating %3.3f GB memory ... ", (double) (TotVertices+2*numEdges)*8/(1<<30));
#endif
	permV = (LONG_T *) malloc(TotVertices*sizeof(LONG_T));
	startVertex = (LONG_T *) malloc(numEdges*sizeof(LONG_T));
	endVertex = (LONG_T *) malloc(numEdges*sizeof(LONG_T));
#if DEBUG
	fprintf(stderr, "done] ");
#endif
	for(i=0; i<TotVertices; i++) {
		permV[i] = i;
	}

	/* Permute the vertices and store them in permV */
	for (i=0; i<TotVertices; i++) {

		t1 = i + sprng(stream) * (TotVertices - i);

		if (t1 != i) {
			t2 = permV[t1];
			permV[t1]=permV[i];
			permV[i]=t2;
		}
	
	}

	for (i=0; i<numEdges; i++) {
		startVertex[i] = permV[startV[i]];
		endVertex[i] = permV[endV[i]];
	}

	free(startV);
	free(endV);	
	free(permV);

	fprintf(stderr, "done\n");

	/*----------------------------------------------*/
        /*              generate edge weights           */
        /*----------------------------------------------*/

        fprintf(stderr, "Generating edge weights ... ");

        weights = (WEIGHT_T *) malloc(numEdges*sizeof(WEIGHT_T));

        for (i=0; i<numEdges; i++) {
                weights[i] = MinWeight + (WEIGHT_T) (MaxWeight - MinWeight) \
                        *sprng(stream);
        }

        fprintf(stderr, "done\n");
	free(stream);
	
	/*-------------------------------------------------------*/
	/*		sort the edge lists with start		 */	
	/*		vertex as primary key 			 */
	/*---------------------------------------------- --------*/

	if (SORT_EDGELISTS) {

		fprintf(stderr, "Sorting edge list by start vertex ... ");

		if (SORT_TYPE == 0) {
			/* Counting sort */
			countingSort(startVertex, endVertex, weights, numEdges);
		} else {
			/* Heap sort */
			heapSort(startVertex, endVertex, weights, numEdges);
		}
	}
	
	fprintf(stderr, "done\n");
	g->start = startVertex;
	g->end = endVertex;
	g->w = weights;
	g->n = TotVertices;
	g->m = numEdges;	

	/*--------------------------------------*/
	/*		update the log		*/
	/*--------------------------------------*/
	
	outfp = fopen(LOGFILE, "a");
	fprintf(outfp, "Graph statistics\n");
	fprintf(outfp, "****************\n");
	fprintf(outfp, "No. of vertices               : %ld\n", g->n);
	fprintf(outfp, "No. of edges                  : %ld\n", g->m);
	fprintf(outfp, "No. of clusters               : %ld\n", totClusters);
	fprintf(outfp, "No. of intra-cluster edges    : %ld\n", numIntraClusterEdges);
	fprintf(outfp, "No. of inter-cluster edges    : %ld\n", numInterClusterEdges);
	
	fclose(outfp);

}
