#!/bin/bash
#
# Outputs alignment results from DALIGNER (with supplementary information from DAZZ_DB) in bella format,
# ASSUMING alignment results have been generated and the below requirements have been met. 
#
# Requires 
# (1) DBNAME to be set externally to the name of the database with which daligner was run.
# (2) LAS to be set externally to e.g. LAS=$DBNAME".las"
# (3) PYPARSE and PYZIP to specify the absolute paths to parserLAdump.py and zip_tag.py, respectively, (from dibella/bella repo) 
#
# Requires in PATH:
# - python (preferrably >= 3.6),
# - LAmerge, LAdump (from DALIGNER), 
# - DBshow (from DAZZ_DB),
#
# An example setup script can be found in /global/homes/m/mme/cori/run-daligner/translation-environment.sh
#   also provided in this directory as example-translation-environment.sh
# If a script is provided as argument 1, it will be executed in the current environment.
#
# @author mme
#
set -x
if [ -z $1 ]; then
  echo "no environment script provided, proceeding under current settings"
else
  source $1
fi
DUMP=$LAS".dump"
IDS=$LAS".ids"
NAMES=$LAS".names"
BELLA=$LAS".bella"

# WARNING LAmerge only if blocks haven't been merged yet
# WARNING LAmerge might fail merging all the blocks is beyond LAmerge's memory limits
# Example merge: LAmerge $LAS *.las && \
LAdump -cl $DBNAME $LAS > $DUMP  # was using -clo, but -o might be overly restrictive... 
python $PYPARSE $DUMP \
&& DBshow -n $DBNAME $IDS > $NAMES \
&& python $PYZIP $NAMES $IDS $BELLA \
&& mv $BELLA $BELLA".txt" # final output if all commands succeed
