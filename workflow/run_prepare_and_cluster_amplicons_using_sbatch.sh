#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --output=%x.o%j
#SBATCH --error=%x.o%j
#SBATCH --partition your_partition
#SBATCH --account=your_account
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

this_script=$(basename $0)

run_script="run_prepare_and_cluster_amplicons.sh"

bash $run_script
