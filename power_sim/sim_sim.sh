#!/bin/bash
#SBATCH --output=test.out
#SBATCH --job-name=test_power
#SBATCH --partition=accel
#SBATCH --time=30-00:00:00
module purge
module load ~/R/
source ~/setup.sh %
Rscript  ~/do_sim.R