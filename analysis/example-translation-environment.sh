#!/bin/bash
module unload python 
module load python/3.6-anaconda-4.4 # Cori
module list #documentation
set -x

DALIGN=$CSCRATCH/DALIGNER/bin
DAZZDB=$CSCRATCH/DAZZ_DB/bin

PYPARSE=$CSCRATCH/dibella-scripts/parserLAdump.py
PYZIP=$CSCRATCH/dibella-scripts/zip_tag.py

export PATH=$DALIGN:$DAZZDB:$PATH

# don't forget to specify DBNAME for the run!
