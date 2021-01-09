COMPILER = g++
CC = gcc
CFLAGS = -I. -O3 -Wall -Wextra -pedantic -ansi -c

optlist.o:	../optlist/optlist.c ../optlist/optlist.h
	$(CC) $(CFLAGS) $<

# flags defined in mtspgemm2017/GTgraph/Makefile.var
check: evaluation.cpp optlist.o 
	$(COMPILER) -o evaluation optlist.o evaluation.cpp

# flags defined in mtspgemm2017/GTgraph/Makefile.var
label: mlEvaluation.cpp optlist.o 
	$(COMPILER) -o mlEvaluation optlist.o mlEvaluation.cpp

clean:
	rm -f *.o
	rm -f check
	rm -f label