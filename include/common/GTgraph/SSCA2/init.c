#include "init.h"

void parseUserInput(int argc, char** argv) {
	int c;
	char configfile[30];
	int outfileSpecified=0;
	int configfileSpecified=0;
	int scaleSpecified=0;	
	if (argc == 1) {
		fprintf(stderr, "GTgraph-ssca2 [-options]\n");
		fprintf(stderr, "\t-s ###  SCALE value (integer) to use (default -- 20)\n"); 
		fprintf(stderr, "\t-c ###  config file to use\n");
		fprintf(stderr, "\t-o ###  output file to write the graph to (default -- sample.gr)\n");
		fprintf(stderr, "\t-h      display this message\n");
		fprintf(stderr, "No config file specified\n");
		fprintf(stderr, "Assigning default values from init.c\n");
		getParams();
		updateLog();

	} else if (argc <= 5) {

		while((c = getopt(argc, argv, "s:c:h:o:")) != -1) {

			switch (c) {
				
			case 's':
				if ((!configfileSpecified) && (!outfileSpecified)) {
					getParams();
				} else if (configfileSpecified) {
					fprintf(stderr, "Warning: Ignoring the SCALE value given as an input argument. The parameters specified in the config file have been applied\n");
					break;
				}
				SCALE = atol(optarg);
                                fprintf(stderr, "Applying a SCALE value of %ld\n", SCALE);
                                TotVertices = 1 << SCALE;
                                MaxCliqueSize = (LONG_T) floor(pow(2, (double) (SCALE/3.0)));
                                MaxWeight = 1 << SCALE;
                                updateLog();
				scaleSpecified = 1;
				break;
					
				case 'o':
				outfileSpecified = 1;
				if ((!configfileSpecified) && (!scaleSpecified)) {
					fprintf(stderr, "No config file specified, assigning default values ...\n");
					getParams();
				}
				WRITE_TO_FILE = 1;
				strcpy(OUTFILE, optarg);
				fprintf(stderr, "Graph will be written to %s\n", OUTFILE);
				updateLog();
				
				break;
				
				case 'c':
				if (scaleSpecified) {
					fprintf(stderr, "Warning: The parameters specified in the config file will be applied and the previous SCALE value will be lost\n");
				}
				configfileSpecified = 1;
				if (!outfileSpecified) {
					if (!scaleSpecified)
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

        fprintf(stderr, "GTgraph-ssca2 [-options]\n");
        fprintf(stderr, "\t-s ###  SCALE value (integer) to use (default -- 20)\n");
        fprintf(stderr, "\t-c ###  config file to use\n");
        fprintf(stderr, "\t-o ###  output file to write the graph to (default -- sample.gr\n");
        fprintf(stderr, "\t-h      display this message\n");

	exit(-1);
}

/* Default Input parameters for graph generation. These values can 
 * also be specified in a configuration file and passed as input to 
 * the graph generator */
void getParams() {

	SCALE = 20;
	TotVertices = (1<<SCALE);
        MaxCliqueSize = (LONG_T) floor(pow(2, (double) (SCALE/3.0)));
	MaxWeight = (1<<SCALE);
	MinWeight = 0;
	
	ProbUnidirectional = 0.2;
	ProbIntercliqueEdges = 0.5;   
	MaxParallelEdges = 3;

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
		else if (strcmp(var,"SCALE")==0) {
			SCALE = (LONG_T) val;
			TotVertices = (1<<SCALE);
		       	MaxCliqueSize = (LONG_T) floor(pow(2, (double) (SCALE/3.0)));
		        MaxWeight = (1<<SCALE);
		} else if (strcmp(var,"TotVertices")==0) {
			TotVertices = (LONG_T) val;
		} else if (strcmp(var,"MaxCliqueSize")==0) {
			MaxCliqueSize = (LONG_T) val;
		} else if (strcmp(var,"MaxWeight")==0) {
			MaxWeight = (WEIGHT_T) val;
		} else if (strcmp(var,"MinWeight")==0) {
			MinWeight = (WEIGHT_T) val;
		} else if (strcmp(var,"ProbUnidirectional")==0) {
			ProbUnidirectional = (double) val;
		} else if (strcmp(var,"ProbInterClEdges")==0) {
			ProbIntercliqueEdges = (double) val;
		} else if (strcmp(var,"MaxParallelEdges")==0) {
			MaxParallelEdges = (LONG_T) val;
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
