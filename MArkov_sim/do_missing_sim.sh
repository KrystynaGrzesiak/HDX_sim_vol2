#!/bin/bash
#SBATCH --output=missing.out
#SBATCH --job-name="missing_sim"
#SBATCH --partition=accel
#
#SBATCH --time=30-00:00:00

module purge

module load system/R−3.6.1

~/R-3.6.1/bin/Rscript simulate_missing_using_markov.R
