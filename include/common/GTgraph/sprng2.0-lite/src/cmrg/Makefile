include ../../../Makefile.var

SRCDIR = ..

all : cmrg.o

cmrg.o : $(SRCDIR)/interface.h cmrg.c  cmrg.h $(SRCDIR)/memory.h  \
		$(SRCDIR)/primes_64.h  $(SRCDIR)/multiply.h \
		$(SRCDIR)/store.h $(SRCDIR)/fwrap_.h
	$(CC) -c $(CFLAGS) $(FFXN) $(INLINEOPT)  cmrg.c -I$(SRCDIR)

clean :
	rm -f *.o *.i

