# README #

Repository containing the source code of the k-mers analysis for long read assembly.

To compile: g++ -std=c++0x trie-tree-implementation.cpp -o trie
To run: ./trie mer_counts_dumps.fa (generated using jellyfish dump) my_genome.fna kmer_length (e.g. 15)

To compile occurrence-matrix.cpp: make parse
To run: ./parse dataset15.fa my_genome.fna 15 (k-mer length) listofdata.txt