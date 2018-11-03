#!/bin/bash
module unload python 
module load python/3.6-anaconda-4.4 # Cori
module list #documentation
set -x

# path should contain path to DALIGNER and DAZZ_DB binaries
SOFT=/global/common/software/m2865/bella-proj/bin
export PATH=$SOFT:$PATH

PYPARSE=$HOME/bella/analysis/parserLAdump.py
PYZIP=$HOME/bella/analysis/zip_tag.py

# don't forget to specify DBNAME for the run!

