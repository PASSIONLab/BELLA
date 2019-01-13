RMATPATH = mtspgemm2017/GTgraph/R-MAT
SPRNPATH = mtspgemm2017/GTgraph/sprng2.0-lite
SEQANPATH = seqan
include mtspgemm2017/GTgraph/Makefile.var
INCLUDE = -I$(SPRNPATH)/include #-I$(CURDIR)/libgaba
SEQINCLUDE = -I$(SEQANPATH)
MLKINCLUDE = -I/opt/intel/composer_xe_2015.0.039/mkl/include
LIBPATH = -L/opt/intel/composer_xe_2015.0.039/mkl/lib 
COMPILER = g++
CC = gcc
CFLAGS = -I. -O3 -W -Wall -Wextra -pedantic -ansi -c
#SEQFLAGS = -DSEQAN_ARCH_SSE4=1 -DSEQAN_ARCH_AVX2=1 -DSEQAN_BGZF_NUM_THREADS=1

sprng:	
	(cd $(SPRNPATH); $(MAKE); cd ../..)

rmat:	sprng
	(cd $(RMATPATH); $(MAKE); cd ../..)

LIBS = -L$(CURDIR)/libbloom/build -lbloom #-L$(CURDIR)/libgaba -lgaba

Buffer.o: kmercode/Buffer.c
	$(CC) -O3 -fopenmp -c -o Buffer.o kmercode/Buffer.c

bound.o: kmercode/bound.cpp
	$(CC) -O3 -fopenmp -c -o bound.o kmercode/bound.cpp

fq_reader.o: kmercode/fq_reader.c
	$(CC) -O3 -std=gnu99 -fopenmp -c -o fq_reader.o kmercode/fq_reader.c

hash_funcs.o: kmercode/hash_funcs.c
	$(CC) -O3 -fopenmp -c -o hash_funcs.o kmercode/hash_funcs.c

bloomlib:
	$(MAKE) -C libbloom all

#gabalib:
#	$(MAKE) -C libgaba all

optlist.o:	optlist/optlist.c optlist/optlist.h
	$(CC) $(CFLAGS) $<

Kmer.o:	kmercode/Kmer.cpp
	$(COMPILER) -fopenmp -std=c++11 -O3 -c -o Kmer.o kmercode/Kmer.cpp

# flags defined in mtspgemm2017/GTgraph/Makefile.var
bella: main.cpp hash_funcs.o fq_reader.o Buffer.o Kmer.o bound.o optlist.o rmat bloomlib 
	#gabalib
	$(COMPILER) -std=c++14 -O3 $(INCLUDE) -march=native -fopenmp -fpermissive $(SEQINCLUDE) -o bella hash_funcs.o Kmer.o Buffer.o fq_reader.o bound.o optlist.o main.cpp ${LIBS}
# add -D__LIBCUCKOO_SERIAL to run lubcuckoo in a single thread
clean:
	(cd mtspgemm2017/GTgraph; make clean; cd ../..)
	rm -f *.o
	rm -f bella
	$(MAKE) -C libbloom clean
	$(MAKE) -C libgaba clean
