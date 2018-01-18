# README #

Repository containing BELLA source code, a long read to long read aligner.

The master currently contains the multithreaded code, while the branch "serial-bella" owns the serial code.

To compile BELLA: make parse
To run BELLA: ./parse <kmers-file> <listoffastq.txt>

In multiply.h: if multitheaded version is used, remember to define which OS you're using (OSX or LINUX) in order to get the memory usage.

The repository contains also the code to get the sensitivy of BELLA as well as other long-read aligners.

To compile:  g++ mtspgemm2017/benchmark.cpp -o <exe>
Tu run: ./exe <ground-truth-file> <bella-output> [<minimap-output>] [<mhap-output>] [<blasr-output>] [<daligner-output>]