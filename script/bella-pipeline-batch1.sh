#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -J bella-pipeline-batch1
#SBATCH --mail-user=gguidi@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -t 30:00

#OpenMP settings:
export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:

# <ksize> <window> <lower> <upper> <xdrop> 
# this is the regular k-mer version
# srun -n 1 -c 64 --cpu_bind=cores ./run-bella-pipeline.sh 17 0 2 7 25

for k in 15 17
do
	for w in 1 2 3 4 5 6
	do
		# <ksize> <window> <lower> <upper> <xdrop> 
		srun -n 1 -c 64 --cpu_bind=cores ./run-bella-pipeline.sh $k $w 2 7 25
	done
done
