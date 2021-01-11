include ../../../Makefile.var

SRCDIR = ..

all : lcg.o 

lcg.o : $(SRCDIR)/interface.h lcg.c  lcg.h $(SRCDIR)/memory.h  \
		$(SRCDIR)/primes_32.h  $(SRCDIR)/multiply.h \
		$(SRCDIR)/store.h $(SRCDIR)/fwrap_.h
	$(CC) -c $(CFLAGS)  $(FFXN) $(INLINEOPT) lcg.c -I$(SRCDIR)

clean :
	rm -f *.o *.f *~ *.i core a.out

