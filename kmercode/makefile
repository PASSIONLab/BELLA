Buffer.o: Buffer.c
	gcc-6 -O3 -c -o Buffer.o Buffer.c

fq_reader.o: fq_reader.c
	gcc-6 -O3 -c -o fq_reader.o fq_reader.c

hash_funcs.o: hash_funcs.c
	gcc-6 -O3 -c -o hash_funcs.o hash_funcs.c

Kmer.o:	Kmer.cpp
	g++-6 -std=c++11 -O3 -c -o Kmer.o Kmer.cpp

parse: main.cpp hash_funcs.o fq_reader.o Buffer.o Kmer.o
	g++-6 -std=c++11 -O3 -o parse hash_funcs.o Kmer.o Buffer.o fq_reader.o main.cpp

clean:
	rm -f *.o
	rm -f parse
