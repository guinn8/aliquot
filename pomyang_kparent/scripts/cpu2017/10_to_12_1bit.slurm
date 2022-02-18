#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gavin.guinn@ucalgary.ca
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --time=05:00:00
#SBATCH --mem=160G
#SBATCH --partition=cpu2017-bf05

../../cli --bound=$((10**12)) --seg_len=$((5 * $((10**8)))) --num_locks=$((10**9)) --preimage_count_bits=1 --num_thread=14
