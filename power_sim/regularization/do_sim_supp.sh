#!/bin/bash
#SBATCH --output=supplementary_pep.out
#SBATCH --job-name=supp
#SBATCH --partition=accel
#SBATCH --time=30-00:00:00
module purge
module load system/R−3.6.1

~/R-3.6.1/bin/Rscript do_sim_supp.R