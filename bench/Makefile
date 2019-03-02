COMPILER = g++
CC = gcc
CFLAGS = -I. -O3 -Wall -Wextra -pedantic -ansi -c 
OMPFLAG = -fopenmp -fpermissive

optlist.o: ../optlist/optlist.c ../optlist/optlist.h
	$(CC) $(CFLAGS) $<

# flags defined in mtspgemm2017/GTgraph/Makefile.var
bench: bench.cpp optlist.o 
	$(COMPILER) -o bench optlist.o bench.cpp

# flags defined in mtspgemm2017/GTgraph/Makefile.var
paf: lostintranslation.cpp optlist.o 
	$(COMPILER) $(OMPFLAG) -o paf optlist.o lostintranslation.cpp

clean:
	rm -f *.o
	rm -f bench
	rm -f paf
