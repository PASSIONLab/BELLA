RMATPATH = GTgraph/R-MAT
SPRNPATH = GTgraph/sprng2.0-lite
include GTgraph/Makefile.var
INCLUDE = -I$(SPRNPATH)/include 
MLKINCLUDE = -I/opt/intel/composer_xe_2015.0.039/mkl/include
LIBPATH = -L/opt/intel/composer_xe_2015.0.039/mkl/lib
COMPILER = g++-4.9

sprng:	
	(cd $(SPRNPATH); $(MAKE); cd ../..)

rmat:	sprng
	(cd $(RMATPATH); $(MAKE); cd ../..)

TOCOMPILE = $(RMATPATH)/graph.o $(RMATPATH)/utils.o $(RMATPATH)/init.o $(RMATPATH)/globals.o 

# flags defined in GTgraph/Makefile.var
mlk_free_spgemm: mkl_free_main.cpp CSC.h CSC.cpp CSR.h CSR.cpp intervaltree/IntervalTree.h multiply.h rmat true-overlap.h compute-delta.h 
	$(COMPILER) $(FLAGS) $(INCLUDE) -fopenmp -o spgemm mkl_free_main.cpp  ${TOCOMPILE} ${LIBS}

clean:
	(cd GTgraph; make clean; cd ../..)
	rm -rf mlk_free_spgemm
