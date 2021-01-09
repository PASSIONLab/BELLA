#!/bin/bash -l
#
# Example slurm script for running daligner with celegans (multiple input files)
# on Cori Haswell.
#
#SBATCH -N 1
#SBATCH -t 04:00:00
#SBATCH --time-min 00:00:30
#SBATCH -L SCRATCH
#SBATCH -C haswell

set -x
LOGICAL_CORES=32
DATA_PATH=$CSCRATCH/pacbio_celegans
DATA_FILES=${DATA_PATH}/all_fastq.txt
DATA=celegans

K="" #"-k17" # form is e.g. -k17
L="" # form is e.g. -l2000  # local alignments will involve 2000 base pairs
DALIGNER_PARAMS="${K} ${L} -A -T"${LOGICAL_CORES}  
RUNDIR=$CSCRATCH/daligner${K}${L}_j${SLURM_JOB_ID}

# should update PATH to include bins for DAZZ_DB, DALIGNER, and seqtk
ENV_SCRIPT=~/bella/analysis/daligner_scripts/example-translation-environment.sh
source $ENV_SCRIPT

mkdir $RUNDIR
cd ${RUNDIR}
cp $0 ./ #documentation

#only need to run the dazz_db setup once per directory
cat ${DATA_FILES} |
while read NAME; do
  seqtk seq -a ${DATA_PATH}/${NAME} > ${NAME}.fasta # changed to -A in seqtk v1.2-8-gd210c57
  fasta2DB ${DATA} ${NAME}.fasta
done
DBsplit ${DATA}

HPC.daligner ${DALIGNER_PARAMS} ${DATA} > align.sh && \
time sh align.sh

