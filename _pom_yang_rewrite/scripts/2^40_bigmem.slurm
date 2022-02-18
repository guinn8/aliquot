#!/bin/bash
#SBATCH --mail-user=gavin.guinn@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --mem=2000G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH --partition=bigmem
#SBATCH --time=48:00:00
\time ../main_prod --bound=$((2**40)) --seg_len=$((2**27)) --writebuf_len=50000000 --preimage_count_bits=8 --num_threads=80

