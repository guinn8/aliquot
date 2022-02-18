#!/bin/bash
#SBATCH --mail-user=gavin.guinn@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --mem=180G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=cpu2021-bf24
#SBATCH --time=24:00:00
\time ../main_prod --bound=$((2**40)) --seg_len=$((2**26)) --writebuf_len=50000000 --preimage_count_bits=1 --num_threads=48

