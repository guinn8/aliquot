#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gavin.guinn@ucalgary.ca
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --time=24:00:00
#SBATCH --mem=800G
#SBATCH --partition=bigmem

# Measured memory util = 544gb
# Estimated util = 565gb
../../cli --bound=$((10**12)) --seg_len=$((1 * $((10**8)))) --num_locks=$((10**9)) --preimage_count_bits=8 --num_thread=72
