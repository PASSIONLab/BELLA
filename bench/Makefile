COMPILER = g++
CC = gcc
CFLAGS = -I. -O3 -Wall -Wextra -pedantic -ansi -c

optlist.o:	../optlist/optlist.c ../optlist/optlist.h
	$(CC) $(CFLAGS) $<

# flags defined in mtspgemm2017/GTgraph/Makefile.var
bench: bench.cpp optlist.o 
	$(COMPILER) -o bench optlist.o bench.cpp

clean:
	rm -f *.o
	rm -f bench
