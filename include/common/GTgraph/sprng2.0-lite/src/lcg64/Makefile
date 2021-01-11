include ../../../Makefile.var
SRCDIR = ..

all : lcg64.o

lcg64.o : $(SRCDIR)/interface.h lcg64.c  lcg64.h $(SRCDIR)/memory.h  \
		$(SRCDIR)/primes_64.h  $(SRCDIR)/multiply.h \
		$(SRCDIR)/store.h $(SRCDIR)/fwrap_.h
	$(CC) -c $(CFLAGS) $(FFXN) $(INLINEOPT) -I$(SRCDIR) lcg64.c 

clean :
	rm -f *.o *.f *~ *.i core a.out
