COMPILER = g++
CC = gcc 
CFLAGS = -I. -O3 -Wall -Wextra -pedantic -ansi -c

optlist.o:	../optlist/optlist.c ../optlist/optlist.h
	$(CC) $(CFLAGS) $<

# flags defined in mtspgemm2017/GTgraph/Makefile.var
benchmark: benchmark.cpp optlist.o 
	$(COMPILER) -o benchmark optlist.o benchmark.cpp

clean:
	rm -f *.o
	rm -f benchmark