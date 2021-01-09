#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -J bella-pipeline-batch2
#SBATCH --mail-user=gguidi@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -t 30:00

#OpenMP settings:
export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:

# Winnowmap requires an additional parameter specifying a file
# containing the list of highly-repetitive k-mers for a reference genome.
# In our experiments, this file was assumed to be pre-computed. Given
# the latest advances in k-mer counting algorithms, this is a low-overhead
# O(n) operation, especially for the values of k typically used for long-read
# mapping (e.g. 15 for ONT reads and 19 for PacBio reads).

for k in 17 15
do
	for w in 1 2 3 4 5 6 7 8 9 10 
	do
		# <ksize> <window> <lower> <upper> <xdrop> 
		srun -n 1 -c 64 --cpu_bind=cores ./run-bella-pipeline2.sh $k $w 2 7 25
	done
done
