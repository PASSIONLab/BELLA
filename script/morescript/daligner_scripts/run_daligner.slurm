#!/bin/bash -l
#
# Example slurm script for running daligner on Cori Haswell.
#
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH --time-min 00:30:00
#SBATCH -L SCRATCH
#SBATCH -C haswell
set -x
LOGICAL_CORES=64
DATA_PATH=$CSCRATCH/pacbio_ecoli30x_ncbi #/pacbio_ecoli30x/
DATA=reads-ecoli-30x-full #pacbio_filtered 
K="" #"-k17" # form is e.g. -k17
L="" # form is e.g. -l2000 # local alignments will involve 2000 base pairs
DALIGNER_PARAMS="${K} ${L} -A -T"${LOGICAL_CORES} 
RUNDIR=$CSCRATCH/daligner${K}${L}_j${SLURM_JOB_ID}

# ENV_SCRIPT should update PATH to include seqtk, DAZZ_DB, and DALIGNER binaries
ENV_SCRIPT=~/bella/analysis/daligner_scripts/example-translation-environment.sh
source $ENV_SCRIPT

mkdir $RUNDIR
cd ${RUNDIR}
cp $0 ./ #documentation

#only need to run the dazz_db setup once per directory
seqtk seq -a ${DATA_PATH}/${DATA}.fastq > ${DATA}.fasta 
fasta2DB ${DATA} ${DATA}.fasta
DBsplit ${DATA}

HPC.daligner ${DALIGNER_PARAMS} ${DATA} > align.sh && \
time sh align.sh

# The following can now be run (ideally in another batch job) to extract and reformat daligner's alignments into bella's format 
#DBNAME=${DATA} PATHTRANS=$HOME/bella/analysis/daligner_scripts/translateDalignerOut.sh translation-environment.sh

