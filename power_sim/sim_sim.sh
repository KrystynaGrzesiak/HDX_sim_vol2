#!/bin/bash
#SBATCH --output=test.out
#SBATCH --job-name=test_power
#SBATCH --partition=accel
#SBATCH --time=30-00:00:00
module purge
module load system/R−3.6.1

~/R-3.6.1/bin/Rscript do_sim_powerHDX.R