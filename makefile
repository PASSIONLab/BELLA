RMATPATH = mtspgemm2017/GTgraph/R-MAT
SPRNPATH = mtspgemm2017/GTgraph/sprng2.0-lite
SEQANPATH = seqan-old
include mtspgemm2017/GTgraph/Makefile.var
INCLUDE = -I$(SPRNPATH)/include
SEQINCLUDE = -I$(SEQANPATH)
MLKINCLUDE = -I/opt/intel/composer_xe_2015.0.039/mkl/include
LIBPATH = -L/opt/intel/composer_xe_2015.0.039/mkl/lib
COMPILER = g++

sprng:	
	(cd $(SPRNPATH); $(MAKE); cd ../..)

rmat:	sprng
	(cd $(RMATPATH); $(MAKE); cd ../..)

TOCOMPILE = $(RMATPATH)/graph.o $(RMATPATH)/utils.o $(RMATPATH)/init.o $(RMATPATH)/globals.o 

Buffer.o: kmercode/Buffer.c
	gcc -g -c -o Buffer.o Buffer.c

fq_reader.o: kmercode/fq_reader.c
	gcc -std=gnu99 -g -c -o fq_reader.o fq_reader.c

hash_funcs.o: kmercode/hash_funcs.c
	gcc -g -c -o hash_funcs.o hash_funcs.c

Kmer.o:	kmercode/Kmer.cpp
	$(COMPILER) -std=c++11 -g -c -o Kmer.o Kmer.cpp

# flags defined in mtspgemm2017/GTgraph/Makefile.var
bella: main.cpp hash_funcs.o fq_reader.o Buffer.o Kmer.o rmat
	$(COMPILER) -std=c++14 $(INCLUDE) -g -O3 -fopenmp -fpermissive $(SEQINCLUDE) -o bella hash_funcs.o Kmer.o Buffer.o fq_reader.o main.cpp ${TOCOMPILE} ${LIBS}

clean:
	(cd mtspgemm2017/GTgraph; make clean; cd ../..)
	rm out.bella
	rm -f *.o
	rm -f bella
