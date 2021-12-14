#!/bin/bash
#SBATCH --output=sim_ibb_models.out
#SBATCH --job-name=sim_ibb_models
#SBATCH --partition=accel
#SBATCH --time=30-00:00:00
module purge
module load system/R−3.6.1

~/R-3.6.1/bin/Rscript do_sim.R