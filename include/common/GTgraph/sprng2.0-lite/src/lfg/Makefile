include ../../../Makefile.var
SRCDIR = ..

all : lfg.o 

lfg.o : $(SRCDIR)/interface.h lfg.c  lfg.h $(SRCDIR)/multiply.h  $(SRCDIR)/memory.h $(SRCDIR)/store.h $(SRCDIR)/fwrap_.h
	$(CC) -c $(CFLAGS)  $(FFXN) $(INLINEOPT) lfg.c -I$(SRCDIR)

clean :
	rm -f *.o *.i

realclean :
	rm -f *.o *.f *~ *.i core a.out

