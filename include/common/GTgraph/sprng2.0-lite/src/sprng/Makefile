include ../../../Makefile.var
SRCDIR = ..
GENLIST = lfg lcg lcg64 cmrg mlfg 

all : sprng 

#---------------------------------------------------------------------------
# Construct Library
#---------------------------------------------------------------------------
SPRNG_COMMON_DEPEND = $(SRCDIR)/primes_32.o $(SRCDIR)/primes_64.o \
			  $(SRCDIR)/fwrap_mpi.o  $(SRCDIR)/cputime.o \
		      $(SRCDIR)/makeseed.o  $(SRCDIR)/simple_mpi.o \
		      $(SRCDIR)/memory.o  $(SRCDIR)/communicate.o  \
			  $(SRCDIR)/checkid.o $(SRCDIR)/store.o

sprng :
	(cd ..; $(MAKE) sprng_common)
	@for g in $(GENLIST) ; do \
		(cd ../$$g; $(MAKE) ) ; \
	done
	rm -f ../../lib/libsprng.a
	$(MAKE) libsprng

libsprng : sprng.o   \
		$(SRCDIR)/lcg/lcg.o $(SRCDIR)/lfg/lfg.o \
		$(SRCDIR)/lcg64/lcg64.o \
		$(SRCDIR)/mlfg/mlfg.o $(SRCDIR)/cmrg/cmrg.o \
		$(SPRNG_COMMON_DEPEND)
	$(AR) $(ARFLAGS) ../../lib/libsprng.a sprng.o \
		$(SRCDIR)/lcg/lcg.o $(SRCDIR)/lfg/lfg.o \
		$(SRCDIR)/lcg64/lcg64.o \
		$(SRCDIR)/mlfg/mlfg.o $(SRCDIR)/cmrg/cmrg.o \
		$(SPRNG_COMMON_DEPEND)
	$(RANLIB) ../../lib/libsprng.a

sprng.o : $(SRCDIR)/interface.h sprng.c   $(SRCDIR)/memory.h  \
	$(SRCDIR)/multiply.h $(SRCDIR)/store.h $(SRCDIR)/fwrap_.h
	$(CC) -c $(CFLAGS)  $(FFXN) $(INLINEOPT) sprng.c -I$(SRCDIR)


.F.f : 
	@if [ -f $*.i ] ; then \
	rm $*.i ;\
	fi
	$(CPP) $*.F
	@if [ -f $*.i ] ; then \
	mv $*.i  $*.f ;\
	fi

#---------------------------------------------------------------------------
clean :
	rm -f *.o *.i

realclean :
	rm -f *.o *.f *~ *.i core a.out

.SUFFIXES : 

.SUFFIXES : .f .F
