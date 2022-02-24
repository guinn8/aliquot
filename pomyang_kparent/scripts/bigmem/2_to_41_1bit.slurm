#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gavin.guinn@ucalgary.ca
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --time=24:00:00
#SBATCH --mem=400G
#SBATCH --partition=bigmem

../../cli --bound=$((2**41)) --seg_len=$((1 * $((2**28)))) --num_locks=$((10**9)) --preimage_count_bits=1 --num_thread=72
