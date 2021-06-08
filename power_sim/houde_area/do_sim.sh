#!/bin/bash
#SBATCH --output=houde_sroude.out
#SBATCH --job-houde_sroude=lme
#SBATCH --partition=accel
#SBATCH --time=30-00:00:00
module purge
module load system/R−3.6.1

~/R-3.6.1/bin/Rscript do_sim.R