#!/bin/bash
#SBATCH --mem=1000000
#SBATCH --nodes=1  
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=75
#SBATCH --partition=bigmem
#SBATCH --time=024:00:00
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./py $1

