#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gavin.guinn@ucalgary.ca
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --time=24:00:00
#SBATCH --mem=700G
#SBATCH --partition=bigmem

../../cli --bound=$((2**42)) --seg_len=$((1 * $((2**29)))) --num_locks=$((10**9)) --preimage_count_bits=1 --num_thread=72
