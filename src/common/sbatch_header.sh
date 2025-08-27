#! /bin/bash
#SBATCH -A your_allocation
#SBATCH -p genomics
#SBATCH --job-name="halla"
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com
#SBATCH --output="halla/logs/halla-%j.out"

#cd $SLURM_SUBMIT_DIR
# 24 hrs, 80Gb and 10 cpus gets job done for most quickly
# Load Conda Environment with Snakemake
module purge all
module load mamba
#mamba init
# source activate halla_biobakery
module load halla
