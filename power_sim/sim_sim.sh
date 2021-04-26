#!/bin/bash
#SBATCH --output=power.out
#SBATCH --job-name=power
#SBATCH --partition=accel
#SBATCH --time=30-00:00:00
module purge
module load system/R−3.6.1

~/R-3.6.1/bin/Rscript do_sim_powerHDX.R