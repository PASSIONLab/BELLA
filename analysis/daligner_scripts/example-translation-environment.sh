#!/bin/bash
module unload python 
module avail python
module load python/3.6-anaconda-5.2 #Cori
module list #documentation 
set -x

# path should contain path to DALIGNER and DAZZ_DB binaries
SOFT=/global/common/software/m2865/bella-proj/bin
export PATH=$SOFT:$PATH

REPO=~/bella
PYPARSE=$REPO/analysis/daligner_scripts/parserLAdump.py
PYZIP=$REPO/analysis/daligner_scripts/zip_tag.py

# don't forget to specify DBNAME for the run!

