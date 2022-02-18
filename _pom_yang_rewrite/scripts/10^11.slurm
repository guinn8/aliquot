#!/bin/bash
#SBATCH --mail-user=gavin.guinn@ucalgary.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mem=60G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --partition=cpu2019-bf05
#SBATCH --time=5:00:00
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
\time ../main $((10**11)) $((5 * $((10 ** 7))))

