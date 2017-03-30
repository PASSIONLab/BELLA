#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <fstream>
#include <iterator>
#include <ctime>
#include <cmath>
#include <string>

//#include "overridenew.h"
#include "utility.h"
#include "longCSC.h"
#include "CSR.h"
#include "IO.h"
#include "multiply.h"
#include "bloom_filter.hpp"
using namespace std;


extern "C" {
#include "GTgraph/R-MAT/defs.h"
#include "GTgraph/R-MAT/init.h"
#include "GTgraph/R-MAT/graph.h"
}

#define VALUETYPE double
#define INDEXTYPE int
//#define GEMM_DEBUG
#define ITERS 10

enum generator_type
{
	rmat_graph,
	er_graph,
};



int main(int argc, char* argv[])
{
	bool binary = false;
	bool product = false;
	bool gen = false;
	string inputname1, inputname2, outputname;

	int edgefactor, scale;
	generator_type gtype;

	if(argc < 5)
	{
		cout << "Normal usage: ./spgemm {gen|binary|text} {rmat|er|matrix1.txt} {scale|matrix2.txt} {edgefactor|product.txt|noout} <numthreads>" << endl;
		return -1;
	}
    	else if(argc < 6)
    	{
		cout << "Normal usage: ./spgemm {gen|binary|text} {rmat|er|matrix1.txt} {scale|matrix2.txt} {edgefactor|product.txt|noout} <numthreads>" << endl;
        	cout << "Running on all " << omp_get_max_threads() << " processors" << endl;
    	}
	else
	{
        	cout << "Running on " << argv[5] << " processors" << endl;
        	omp_set_num_threads(atoi(argv[5]));
	}

	if(string(argv[1]) == string("gen"))
	{
		gen = true;
		cout << "Using synthetically generated data, with " << argv[2] << " generator of scale " << argv[3] << " and edgefactor " <<  argv[4] << endl;
		scale = atoi(argv[3]);
		edgefactor = atoi(argv[4]);
		if(string(argv[2]) == string("rmat"))
			gtype = rmat_graph;
		else
			gtype = er_graph;
	}
	else
	{
		inputname1 =  argv[2];
               	inputname2 =  argv[3];
		string isbinary(argv[1]);
               	if(isbinary == "text")
                       	binary = false;
             	else if(isbinary == "binary")
                       	binary = true;
               	else
                       	cout << "unrecognized option, assuming text file" << endl;

		if(string(argv[4]) != string("noout"))
		{
			product = true;
			outputname = argv[4];
		}
	}

	CSC<INDEXTYPE,VALUETYPE> * A_csc, * B_csc, * C_csc_verify;
	if(gen)
	{
		double a, b, c, d;
		if(gtype == rmat_graph)
		{
			a = 0.45;
        		b = 0.15;
        		c = 0.15;
        		d = 0.25;
		}
		else
		{
			a = b =  c = d = 0.25;
		}
		getParams();
		setGTgraphParams(scale, edgefactor, a, b, c, d);
		graph G1;
        	graphGen(&G1);
		cerr << "Generator returned" << endl;
		A_csc = new longCSC<INDEXTYPE,VALUETYPE> (G1);	// convert to CSC
		if (STORE_IN_MEMORY) {
                	free(G1.start);
                	free(G1.end);
                	free(G1.w);
        	}

		graph G2;
        	graphGen(&G2);
		cerr << "Generator returned" << endl;
		B_csc = new longCSC<INDEXTYPE,VALUETYPE> (G2);	// convert to CSC

		if (STORE_IN_MEMORY) {
                	free(G2.start);
                	free(G2.end);
                	free(G2.w);
        	}
		
	}
	else
	{
		if(binary)
		{
        		ReadBinary( inputname1, A_csc );
        		ReadBinary( inputname2, B_csc );
			if(product)
				ReadBinary( outputname, C_csc_verify );
		
		}
		else 
		{
			cout << "reading input matrices in text(ascii)... " << endl;
        		ReadASCII( inputname1, A_csc );
        		ReadASCII( inputname2, B_csc );
			if(product)
				ReadASCII( outputname, C_csc_verify );
		}
	}
	
  	A_csc->Sorted();
  	B_csc->Sorted();
    CSC<INDEXTYPE,VALUETYPE> C_csc;
    
    double start = omp_get_wtime( );

    for(int i=0; i< ITERS; ++i)
        HeapSpGEMM_gmalloc(*A_csc, *B_csc, multiplies<VALUETYPE>(), plus<VALUETYPE>(), myidentity<VALUETYPE>(), C_csc);
    
    cout << "HeapSpGEMM returned with " << C_csc.nnz << " nonzeros" << endl;
    double end = omp_get_wtime( );
    printf("HeapSpGEMM, start = %.16g,nend = %.16g, diff = %.16g\n", start, end, (end - start)/ITERS);

    C_csc.Sorted();
    
    if(product)
	{
		if(*C_csc_verify == C_csc)
			cout << "HeapSpGEMM is correct" << endl;
		else
			cout << "HeapSpGEMM is INcorrect" << endl;

	}
}
