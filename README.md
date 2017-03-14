# README #

Repository containing the source code of the k-mers analysis for long read assembly.

To compile: g++ -std=c++0x trie-tree-implementation.cpp -o trie
To run: ./trie mer_counts_dumps.fa (generated using jellyfish dump) my_genome.fna kmer_length (e.g. 15)
