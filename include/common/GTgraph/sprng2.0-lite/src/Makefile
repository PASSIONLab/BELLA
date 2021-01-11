include ../../Makefile.var

LIBLIST = lfg lcg lcg64 cmrg mlfg 

#---------------------------------------------------------------------------
# Jump point to sprng
#---------------------------------------------------------------------------
all : 
	(cd sprng; $(MAKE) )

#---------------------------------------------------------------------------
# Create files that all generators need
#---------------------------------------------------------------------------
SPRNG_COMMON_DEPEND = primes_32.o primes_64.o \
			fwrap_mpi.o  cputime.o makeseed.o  store.o\
                    simple_mpi.o memory.o  communicate.o  checkid.o

sprng_common : $(SPRNG_COMMON_DEPEND)

primes_32.o : primes_32.h primes_32.c
	$(CC) -c $(CFLAGS) primes_32.c 

primes_64.o : primes_64.h primes_64.c
	$(CC) -c $(CFLAGS) primes_64.c 

fwrap_mpi.o : fwrap.h fwrap_mpi.c
	$(CC) -c $(CFLAGS) $(FFXN) fwrap_mpi.c

cputime.o : cputime.c fwrap.h
	$(CC) -c $(CFLAGS) $(FFXN) cputime.c

makeseed.o : makeseed.c
	$(CC) -c $(CFLAGS) makeseed.c

simple_mpi.o : simple_mpi.c memory.h interface.h
	$(CC) -c $(CFLAGS) simple_mpi.c

memory.o : memory.c
	$(CC) -c $(CFLAGS) memory.c

communicate.o : communicate.c memory.h
	$(CC) -c $(CFLAGS) communicate.c

checkid.o : checkid.c memory.h
	$(CC) -c $(CFLAGS) checkid.c

store.o : store.c store.h
	$(CC) -c $(CFLAGS) store.c

.SUFFIXES : 

.SUFFIXES : .f .F

#---------------------------------------------------------------------------
.F.f : 
	@if [ -f $*.i ] ; then \
	rm $*.i ;\
	fi
	$(CPP) $(SRCDIR) $*.F
	@if [ -f $*.i ] ; then \
	mv $*.i  $*.f ;\
	fi

#---------------------------------------------------------------------------
clean :
	rm -f *.o *.i *~ core a.out
	@for l in $(LIBLIST) sprng;  do \
          cd $$l ; \
          $(MAKE) PLAT=$(PLAT) clean ; \
          cd .. ; \
        done

