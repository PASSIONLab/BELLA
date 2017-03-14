#include "init.h"

void parseUserInput(int argc, char** argv) {
	int c;
	char configfile[30];
	int outfileSpecified = 0;
	int configfileSpecified = 0;
	int typeSpecified = 0;
	if (argc == 1) {
		fprintf(stderr, "GTgraph-random [-options]\n");
		fprintf(stderr, "\t-c ###  config file to use\n");
		fprintf(stderr, "\t-t ###  random graph model\n");
		fprintf(stderr, "\t        (0 for G(n, p), 1 for G(n, m))\n"); 
		fprintf(stderr, "\t-n ###  no. of vertices\n");
		fprintf(stderr, "\t-p ###  probability of an edge (graph model 0)\n");        		 
		fprintf(stderr, "\t-m ###  no. of edges (graph model 1)\n");	
		fprintf(stderr, "\t-o ###  output file to write the graph to\n");
		fprintf(stderr, "\t        (default: sample.gr)\n"); 
		fprintf(stderr, "\t-h      display this message\n\n");
		fprintf(stderr, "No config file specified\n");
		fprintf(stderr, "Assigning default values from init.c\n");
		getParams();
		updateLog();

	} else if (argc <= 9) {
		
		while ((c = getopt(argc, argv, "c:t:n:m:p:o:h:")) != -1) {

			switch (c) {
				
				case 't':
				if (configfileSpecified) 
					break;
				if (!outfileSpecified) {
					getParams();
				}
				GRAPH_MODEL = atoi(optarg);
				typeSpecified = 1;
                                fprintf(stderr, "Graph model %d chosen\n", GRAPH_MODEL);
                                updateLog();
				break;
					
				case 'o':
				outfileSpecified = 1;
				if ((!configfileSpecified) && (!typeSpecified)) {
					fprintf(stderr, "No config file specified, assigning default values ...\n");
					getParams();
				}
				WRITE_TO_FILE = 1;
				strcpy(OUTFILE, optarg);
				fprintf(stderr, "Graph will be written to %s\n", OUTFILE);
				updateLog();
				break;
				
				case 'c':
				fprintf(stderr, "Warning: The parameters specified in the config file will be applied and all other options will be lost\n");
				configfileSpecified = 1;
				if (!outfileSpecified) {
					getParams();
					strcpy(configfile, optarg);
					fprintf(stderr, "Reading config file %s\n", configfile);
					getParamsFromFile(configfile);
				} else {
					strcpy(configfile, optarg);
					fprintf(stderr, "Updating parameters from config file %s\n", configfile);
					getParamsFromFile(configfile);
				}
				updateLog();
				break;
		
				case 'n':
				if (configfileSpecified)
					break;
				if (!typeSpecified) {
					fprintf(stderr, "Error! First specify graph model using the -t option before this argument\n");
					exit(-1);
				}

				n = atol(optarg);
				fprintf(stderr, "n is set to %ld\n", n);
				updateLog();
				break;

				case 'p':
				if (configfileSpecified)
                                        break;
                                if (!typeSpecified) {
                                        fprintf(stderr, "Error! First specify graph model using the -t option before this argument\n");
                                        exit(-1);
                                }

                                p = atof(optarg);
                                fprintf(stderr, "p is set to %lf\n", p);
                                updateLog();
				break;

				case 'm':
				if (configfileSpecified)
					break;
				if (!typeSpecified) {
                                        fprintf(stderr, "Error! First specify graph model using the -t option before this argument\n");
                                        exit(-1);
                                }
				m = atol(optarg);
                                fprintf(stderr, "m is set to %ld\n", m);
                                updateLog();
                                break;

				case 'h':
				usage();

				default:
				usage();
			}
		}

		
	} else {
		fprintf(stderr, "Invalid input arguments\n");
		usage();
	}
	
}

void usage() {

	fprintf(stderr, "GTgraph-random [-options]\n");
        fprintf(stderr, "\t-c ###  config file to use\n");
        fprintf(stderr, "\t-t ###  random graph model\n");
        fprintf(stderr, "\t        (0 for G(n, p), 1 for G(n, m))\n");
        fprintf(stderr, "\t-n ###  no. of vertices\n");
        fprintf(stderr, "\t-p ###  probability of an edge (graph model 0)\n");
        fprintf(stderr, "\t-m ###  no. of edges (graph model 1)\n");
        fprintf(stderr, "\t-o ###  output file to write the graph to\n");
        fprintf(stderr, "\t        (default: no write)\n");
        fprintf(stderr, "\t-h      display this message\n\n");
	
	exit(-1);
}

/* Default Input parameters for graph generation. These values can 
 * also be specified in a configuration file and passed as input to 
 * the graph generator */
void getParams() {

	GRAPH_MODEL = 1;
	n = 10000000;
	p = 0.000001;
	m = 100000000;
	MAX_WEIGHT = 100; 
	MIN_WEIGHT = 0;

	SELF_LOOPS = 0;
	STORE_IN_MEMORY = 1;
	
	SORT_EDGELISTS = 1;
	SORT_TYPE = 0;
	WRITE_TO_FILE = 1;

	strcpy(OUTFILE, "sample.gr");	
	strcpy(LOGFILE, "log");	
}

void getParamsFromFile(char* configfile) {

	/* read parameters from config file */
	FILE *fp;
	char line[128], var[32];
	double val;

	fp = fopen(configfile,"r");
	if (fp == NULL) {
		fprintf(stderr, "Unable to open config file:%s\n",configfile);
		exit(-1);
	}

	while (fgets(line, sizeof (line), fp) != NULL) {
		sscanf(line,"%s %lf",var, &val);
		if (*var == '#') continue;  /* comment */
		else if (strcmp(var,"GRAPH_MODEL")==0) {
			GRAPH_MODEL = (int) val;
		} else if (strcmp(var,"n")==0) {
			n = (LONG_T) val;
		} else if (strcmp(var,"m")==0) {
			m = (LONG_T) val;
		} else if (strcmp(var, "p")==0) {
			p = (double) val;
		} else if (strcmp(var,"MAX_WEIGHT")==0) {
			MAX_WEIGHT = (WEIGHT_T) val;
		} else if (strcmp(var,"MIN_WEIGHT")==0) {
			MIN_WEIGHT = (WEIGHT_T) val;
		} else if (strcmp(var,"STORE_IN_MEMORY")==0) {
			STORE_IN_MEMORY = (int) val;
		} else if (strcmp(var, "SELF_LOOPS")==0) {
			SELF_LOOPS = (int) val;
		} else if (strcmp(var,"SORT_EDGELISTS")==0) {
			SORT_EDGELISTS = (int) val;
		} else if (strcmp(var,"SORT_TYPE")==0) {
			SORT_TYPE = (int) val;
		} else if (strcmp(var,"WRITE_TO_FILE")==0) {
			WRITE_TO_FILE = (int) val;
		} else {
			fprintf(stderr,"unknown parameter: %s\n",line);
		}
	}

}
