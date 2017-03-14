#include "graph.h"

int *stream1, *stream2, *stream3, *stream4;
void graphGen(graph* g)
{

	LONG_T *startVertex, *endVertex;
	WEIGHT_T *weights;
	LONG_T i, u, v, step;
	LONG_T numEdges;
	WEIGHT_T w;
	FILE* outfp;
	double a0, b0, c0, d0;

	/*----------------------------------------------*/
	/*		initialize SPRNG 		*/
	/*----------------------------------------------*/

	stream1 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
	stream2 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED2, SPRNG_DEFAULT);
	stream3 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED3, SPRNG_DEFAULT);
	stream4 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED4, SPRNG_DEFAULT);

	/*------------------------------------------------------*/
	/*		generate edges as per the		*/		 
	/*		graph model and user options	      	*/
	/*------------------------------------------------------*/

	if ((STORE_IN_MEMORY == 0) && (SORT_EDGELISTS == 0)) {
		fprintf(stderr, "Generating edges 'on the fly' and writing to %s\n", OUTFILE);
		outfp = fopen(OUTFILE, "w");
                fprintf(outfp, "c FILE			: %s\n", OUTFILE);
                fprintf(outfp, "c No. of vertices	: %ld\n", n);
		/* Leave a blank line to fill in the no. of edges later */
		fprintf(outfp, "                                                 \n");
                fprintf(outfp, "c Max. weight	        : %ld\n", MAX_WEIGHT);
                fprintf(outfp, "c Min. weight          	: %ld\n", MIN_WEIGHT);
                fprintf(outfp, "c A directed arc from u to v of weight w\n");
                fprintf(outfp, "c is represented below as ' a  u  v  w '\n");
		fprintf(stderr, "Generating Edges ... ");
		fprintf(outfp, "p sp %ld %ld\n", n, m);

		a0 = a; b0 = b; c0 = c; d0 = d;
		for (i=0; i<m; i++) {
			a = a0; b = b0; c = c0; d = d0;
			u = 1;
			v = 1;
			step = n/2;

			while (step >= 1) {
				choosePartition(&u, &v, step);
				step = step/2;
				varyParams(&a, &b, &c, &d);	
			}

			/* Create edge [u-1, v-1] */
			if ((!SELF_LOOPS) && (u == v)) {
				i--;
				continue;
			}

			w = MIN_WEIGHT + (WEIGHT_T) (MAX_WEIGHT - MIN_WEIGHT) * sprng(stream2);
                        fprintf(outfp, "a %ld %ld %ld\n", u, v, w);
			
		}

		fclose(outfp);
                fprintf(stderr, "done\n");

                free(stream1);
                free(stream2);
		free(stream3);
		free(stream4);
                return;

	}

	
	fprintf(stderr, "Generating edges ... ");

	startVertex = (LONG_T *) malloc(m * sizeof(LONG_T));
	endVertex   = (LONG_T *) malloc(m * sizeof(LONG_T));

	a0 = a; b0 = b; c0 = c; d0 = d;
	for (i=0; i<m; i++) {
		a = a0; b = b0; c = c0; d = d0;
	        u = 1;
                v = 1;
                step = n/2;
		
                while (step >= 1) {
        		choosePartition(&u, &v, step);
                        step = step/2;
                        varyParams(&a, &b, &c, &d);
                }

                /* Create edge [u-1, v-1] */
                if ((!SELF_LOOPS) && (u == v)) {
                        i--;
                        continue;
                }
		
		startVertex[i] = u-1;
		endVertex[i] = v-1;

        }

	numEdges = m;
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
                if (SORT_TYPE == 0) {
                       	/* Counting sort */
                       	countingSort(startVertex, endVertex, weights, numEdges);
                } else {
                       	/* Heap sort */
                   	heapSort(startVertex, endVertex, weights, numEdges);
                }
                fprintf(stderr, "done\n");
        }

        g->start = startVertex;
        g->end = endVertex;
        g->w = weights;
        g->n = n;
        g->m = numEdges;

	free(stream3);
	free(stream4);

}


void choosePartition(LONG_T *u, LONG_T* v, LONG_T step) {

	double p;

	p = sprng(stream1);

	if (p < a) {

		/* Do nothing */

	} else if ((a < p) && (p < a+b)) {

		*v = *v + step;

	} else if ((a+b < p) && (p < a+b+c)) {
	
		*u = *u + step;

	} else if ((a+b+c < p) && (p < a+b+c+d)) {

		*u = *u + step;
		*v = *v + step;
	}

}

void varyParams(double* a, double* b, double* c, double* d) {

	double v, S;

	/* Allow a max. of 5% variation */
	v = 0.05;
	
	if (sprng(stream4) > 0.5)
		*a += *a * v * sprng(stream3);
	else 
		*a -= *a * v * sprng(stream3);

	if (sprng(stream4) > 0.5)
                *b += *b * v * sprng(stream3);
        else 
                *b += *b * v * sprng(stream3);

	if (sprng(stream4) > 0.5)
                *c += *c * v * sprng(stream3);
        else 
                *c -= *c * v * sprng(stream3);

	if (sprng(stream4) > 0.5)
                *d += *d * v * sprng(stream3);
        else 
                *d -= *d * v * sprng(stream3);


	 S = *a + *b + *c + *d;
	*a = *a/S;
	*b = *b/S;
	*c = *c/S;
	*d = *d/S;

}
