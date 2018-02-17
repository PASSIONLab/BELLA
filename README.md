# README #

Repository containing BELLA source code, a long read aligner to de novo genome assembly.

The master contains the multithreaded code on single node.

To compile BELLA: make bella
To run BELLA: ./bella -f <kmers-file> -i <listoffastq> -o <out-filename> [-k <kmer-length>] [-a <alignment-score-thr>] [-p <alignment-xdrop-factor>] [-z]
	-h will show the usage
	
The repository contains also the code to get the sensitivy of BELLA as well as other long-read aligners.

To compile: g++ mtspgemm2017/benchmark.cpp -o <exe>
Tu run: ./exe <ground-truth-file> <bella-output> [<minimap-output>] [<mhap-output>] [<blasr-output>] [<daligner-output>]