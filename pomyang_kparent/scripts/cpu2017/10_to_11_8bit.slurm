#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gavin.guinn@ucalgary.ca
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --time=05:00:00
#SBATCH --mem=100G
#SBATCH --partition=cpu2017-bf05

../../cli --bound=$((10**11)) --seg_len=$((5 * $((10**7)))) --num_locks=$((10**8)) --preimage_count_bits=8 --num_thread=14
