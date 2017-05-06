#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {

	FILE *in;
	FILE *mean;
	FILE *out;
	int depth = atoi(argv[2]);
	int length = atoi(argv[3]);
	int elem, num, i, max=0, RKS=0, TKS=0;						
	char *buf = malloc(sizeof(char)*(length+2)), *res;
	int uRKS = 2*depth, lRKS = depth/2;				
	

	printf("\nThe input file is: %s", argv[1]);
        printf("\nThe k-mer length is: %d", length);
	printf("\nThe pbsim depth is: %d\n", depth);
	
	in = fopen (argv[1], "r");
	if(in == NULL){
		perror("\nFile [in] opening FAILED\n");
		exit(1);
	}

	mean = fopen("first_step.txt", "w+"); //temporary solution to handle the huge amount of data
	if(mean == NULL){
		perror("\nFile [mean] opening FAILED\n");
		exit(1);
	}	

	while(!feof(in)){
		res = fgets(buf, length+2, in);
		if(res == NULL)
			break;
		
		if((buf[0]!='A') && (buf[0]!='C') && (buf[0]!='G') && (buf[0]!='T')) { 		
				sscanf(buf, ">%d", &elem);
				if((elem >= lRKS) && (elem <= uRKS)) { 	
					fprintf(mean, "%d\n", elem);
					RKS++;
					if(elem > max)
						max = elem;
				}
		}
		else
		TKS++;										
	}

	int *count = malloc(sizeof(int)*(max));	

	out = fopen("output.txt", "w");
	if(out == NULL){
		perror("\nFile [out] opening FAILED\n");
		exit(1);
	}

	fseek(mean, 0, SEEK_SET);

	while(!feof(mean)) {
		res = fgets(buf, length+2, mean);
		if(res == NULL)
			break;
				sscanf(buf, "%d", &num);
				count[num]++;	
	}

	for(i=0; i<=max; i++) {
		if(count[i]>0) {
			fprintf(out, "%d %d\n", i, count[i]);
		}
	}

	printf("Total number of k-mers: %d\nReliable k-mers (RKS) number: %d\n\n", TKS, RKS);
	printf("The output.txt file contains the frequencies of the RKS occurrences.\n\n");	
	
	fclose(in);
	fclose(mean);
	fclose(out);
	
	return 0;

}