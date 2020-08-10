#!/bin/bash
#SBATCH --time=7-0:0
#SBATCH --mem=10000              
#SBATCH --nodes=1            
#SBATCH --ntasks=1           
#SBATCH --partition=cpu2019
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=gavinguinn1@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./main
