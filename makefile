RMATPATH = mtspgemm2017/GTgraph/R-MAT
SPRNPATH = mtspgemm2017/GTgraph/sprng2.0-lite
DALPATH = daligner/
include mtspgemm2017/GTgraph/Makefile.var
INCLUDE = -I$(SPRNPATH)/include -I$(DALPATH)
MLKINCLUDE = -I/opt/intel/composer_xe_2015.0.039/mkl/include
LIBPATH = -L/opt/intel/composer_xe_2015.0.039/mkl/lib
COMPILER = g++
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

sprng:	
	(cd $(SPRNPATH); $(MAKE); cd ../..)

rmat:	sprng
	(cd $(RMATPATH); $(MAKE); cd ../..)

align.o: $(DALPATH)/align.c 
	gcc $(CFLAGS) -g -c -o align.o $(DALPATH)/align.c -lpthread -lm

DB.o: $(DALPATH)/DB.c 
	gcc $(CFLAGS) -g -c -o DB.o $(DALPATH)/DB.c -lpthread -lm

QV.o: $(DALPATH)/QV.c 
	gcc $(CFLAGS) -g -c -o QV.o $(DALPATH)/QV.c -lpthread -lm

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
parse: occurrence-matrix.cpp hash_funcs.o fq_reader.o Buffer.o Kmer.o rmat align.o DB.o QV.o
	$(COMPILER) -std=c++11 $(INCLUDE) -O3 -fopenmp -fpermissive -ledlib -o parse hash_funcs.o Kmer.o Buffer.o fq_reader.o align.o DB.o QV.o occurrence-matrix.cpp ${TOCOMPILE} ${LIBS}

clean:
	(cd mtspgemm2017/GTgraph; make clean; cd ../..)
	rm -f *.o
	rm -f parse
