#!/bin/bash
#SBATCH --output=sim_lme.out
#SBATCH --job-name=lme
#SBATCH --partition=accel
#SBATCH --time=30-00:00:00
module purge
module load system/R−3.6.1

~/R-3.6.1/bin/Rscript do_sim.R