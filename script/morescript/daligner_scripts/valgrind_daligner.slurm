#!/bin/bash -l
#
# Example script for Cori Haswell.
# Assumes a directory of daligner results (including align.sh) exists.
# Valgrind will introduce additional memory overhead. 
# 
#SBATCH -N 1
#SBATCH -c 64
#SBATCH -C haswell
#SBATCH --time-min=00:00:30
#SBATCH -t 02:00:00

module load valgrind
set -x

# 
# DIR should be set to the location of existing daligner results (database and align.sh).
# To collect fresh results, see example script run_daligner.slurm.
#
# To run valgrind in a separate directory:
#   mkdir <preferred directory location> \
#   && cd <new directory> \
#   && ln -s $DIR ./
#
DIR=/global/cscratch1/sd/mme/daligner_j16235709 #ecoli 30x (sample)
cd $DIR
cp $0 ./ # copy this script to the run directory (documentation)
#
# PATH should include all tools required by align.sh (most notably daligner compiled with -g for valgrind)
#
export PATH=/global/common/software/m2865/bella-proj/g0-bin:$PATH
#
# generates a new script based on align.sh with valgrind... in place of daligner
#
NEW=valgrind_align.sh
cp align.sh $NEW
sed -i -e 's/daligner/valgrind --tool=massif --pages-as-heap=yes daligner/g' $NEW && cat $NEW
set +x
#
# runs the valgrind script
#
sh $NEW 


