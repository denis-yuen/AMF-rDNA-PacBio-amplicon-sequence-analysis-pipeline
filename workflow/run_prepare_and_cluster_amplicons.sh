#!/bin/bash
# Shell script to run Snakefile with conda

snakefile_name="Snakefile_prepare_and_cluster_amplicons"

script_name=$(basename $0)
today=$(date +"%Y%m%d")
log_file=$script_name.$today."log"

conda_location="/path/to/conda/bin"
export PATH=$conda_location:$PATH

eval "$(conda shell.bash hook)"
conda activate snakemake

## to check DAG before starting run
snakemake --unlock -n -s $snakefile_name all > $log_file 

## to run workflow
snakemake --use-conda --conda-frontend conda --conda-prefix envs -j 4 -s $snakefile_name all


