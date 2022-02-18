#!/bin/bash
#SBATCH --mail-user=gavin.guinn@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --mem=117G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu2013
#SBATCH --time=168:00:00
\time ../main_prod --bound=$((2**40)) --seg_len=$((2**10)) --writebuf_len=100000 --preimage_count_bits=1 --num_threads=16

