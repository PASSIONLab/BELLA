#include "init.h"
#include "math.h"

void parseUserInput(int argc, char** argv) {
	int c;
	char configfile[30];
	int outfileSpecified = 0;
	int configfileSpecified = 0;
	int varSpecified = 0;
	if (argc == 1) {
		fprintf(stderr, "GTgraph-rmat [-options]\n");
		fprintf(stderr, "\t-c ###  config file to use\n");
		fprintf(stderr, "\t-n ###  no. of vertices       (default -- 10^7)\n");
		fprintf(stderr, "\t-m ###  no. of directed edges (default -- 10^8)\n");        		 
		fprintf(stderr, "\t-o ###  output file to write the graph to\n");
		fprintf(stderr, "\t-h      display this message\n\n");
		fprintf(stderr, "No config file specified\n");
		fprintf(stderr, "Assigning default values from init.c\n");
		getParams();
		updateLog();

	} else if (argc <= 9) {

		while((c = getopt(argc, argv, "c:n:m:o:h:")) != -1) {

			switch (c) {
				
				case 'o':
				outfileSpecified = 1;
				if ((!configfileSpecified) && (!varSpecified)) {
					fprintf(stderr, "No config file specified, assigning default values ...\n");
					getParams();
				}
				WRITE_TO_FILE = 1;
				strcpy(OUTFILE, optarg);
				fprintf(stderr, "Graph will be written to %s\n", OUTFILE);
				updateLog();
				break;
				
				case 'c':
				fprintf(stderr, "Warning: The parameters specified in the config file will be applied and all other input arguments will be ignored\n");
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
				if (varSpecified == 0) {
					getParams();
					varSpecified = 1;
				}
				n = atol(optarg);
				fprintf(stderr, "n is set to %ld\n", n);
				updateLog();
				break;

				case 'm':
				if (configfileSpecified)
					break;
				if (varSpecified == 0) {
                                        getParams();
                                        varSpecified = 1;
                                }
                                m = atol(optarg);
                                fprintf(stderr, "m is set to %ld\n", n);
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

	fprintf(stderr, "GTgraph-rmat [-options]\n");
        fprintf(stderr, "\t-c ###  config file to use\n");
        fprintf(stderr, "\t-n ###  no. of vertices       (default -- 10^7)\n");
        fprintf(stderr, "\t-m ###  no. of directed edges (default -- 10^8)\n");
        fprintf(stderr, "\t-o ###  output file to write the graph to\n");
        fprintf(stderr, "\t-h      display this message\n\n");
	
	exit(-1);
}

/* Default Input parameters for graph generation. These values can 
 * also be specified in a configuration file and passed as input to 
 * the graph generator */
void getParams() {

	n = 10000000;
	m = 100000000;

	a = 0.45;
	b = 0.15;
	c = 0.15;
	d = 0.25;

	MAX_WEIGHT = 100; 
	MIN_WEIGHT = 0;

	SELF_LOOPS = 0;
	STORE_IN_MEMORY = 1;
	
	SORT_EDGELISTS = 1;
	SORT_TYPE = 0;
	WRITE_TO_FILE = 0;

	strcpy(OUTFILE, "sample.gr");	
	strcpy(LOGFILE, "log");	
}

void setGTgraphParams(int scale, int edgefactor, double a0, double b0, double c0, double d0)
{
	n = pow(2.0, scale);
	m = n * edgefactor;
	a = a0;
	b = b0;
	c = c0;
	d = d0;
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
		else if (strcmp(var,"n")==0) {
			n = (LONG_T) val;
		} else if (strcmp(var,"m")==0) {
			m = (LONG_T) val;
		} else if (strcmp(var, "a")==0) {
			a = (double) val;
		} else if (strcmp(var, "b")==0) {
			b = (double) val;
		} else if (strcmp(var, "c")==0) {
			c = (double) val;
		} else if (strcmp(var, "d")==0) {
			d = (double) val;
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

	fclose(fp);

}
