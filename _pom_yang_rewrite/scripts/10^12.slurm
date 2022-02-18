#!/bin/bash
#SBATCH --mail-user=gavin.guinn@ucalgary.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mem=1250G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=70
#SBATCH --partition=bigmem
#SBATCH --time=12:00:00
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
\time ../main $((10**12)) $((5 * $((10 ** 8))))

