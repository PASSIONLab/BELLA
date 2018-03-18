# README #

Repository containing BELLA source code, a long read aligner to de novo genome assembly.

The master contains the multithreaded code on single node.
#
To compile BELLA: make bella
#
To run BELLA: ./bella -f <kmers-file> -i <listoffastq> -o <out-filename> [-k <kmer-length>] [-a <alignment-score-thr>] [-p <alignment-xdrop-factor>] [-z]
	-h will show the usage
#	
The repository contains also the code to get the sensitivy of BELLA as well as other long-read aligners.
#
To compile: cd analysis && make check
# 
To run: ./evaluation -g <ground-truth-file> -n <nreads> -b <bella-output> [-m <minimap-output>] [-p <mhap-output>] [-l <blasr-output>] [-d <daligner-output>] [-z]
	-h will show the usage
#
SAMparser.py allow to transform the BWA-MEM .sam outfile in a simpler format usable as input to the evaluation code when using real dataset
#
Requirement: simplesam package, it can be installed via pip: pip install simplesam
To run: python SAMparser.py <bwamem-output>
#
Analysis folder contains also the code related to the Markov model.