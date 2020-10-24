#!/bin/bash
#SBATCH --output=debug.out
#SBATCH --job-name="debug_vec"
#SBATCH --partition=accel
#
#SBATCH --time=30-00:00:00

module purge

module load system/R−3.6.1

~/R-3.6.1/bin/Rscript debug_vec.R
