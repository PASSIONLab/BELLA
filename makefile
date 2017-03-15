Buffer.o: kmercode/Buffer.c
	gcc -O3 -c -o Buffer.o Buffer.c

fq_reader.o: kmercode/fq_reader.c
	gcc -O3 -c -o fq_reader.o fq_reader.c

hash_funcs.o: kmercode/hash_funcs.c
	gcc -O3 -c -o hash_funcs.o hash_funcs.c

Kmer.o:	kmercode/Kmer.cpp
	g++ -std=c++11 -O3 -c -o Kmer.o Kmer.cpp

parse: occurrence-matrix.cpp hash_funcs.o fq_reader.o Buffer.o Kmer.o
	g++ -std=c++11 -O3 -o parse hash_funcs.o Kmer.o Buffer.o fq_reader.o occurrence-matrix.cpp

clean:
	rm -f *.o
	rm -f parse
