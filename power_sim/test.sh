#!/bin/bash
#SBATCH --output=test_power.out
#SBATCH --job-name=test_power_test
#SBATCH --partition=accel
#SBATCH --time=30-00:00:00
module purge
module load system/R−3.6.1

~/R-3.6.1/bin/Rscript test_do_sim_powerHDX.R