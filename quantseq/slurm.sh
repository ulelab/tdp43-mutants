#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=cpu

module purge
module load Anaconda3

source activate tdp43-mutants

mkdir logs

srun snakemake --cluster "sbatch {params.cluster}" --jobs 100
