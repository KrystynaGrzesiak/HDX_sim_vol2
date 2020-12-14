#!/bin/bash
#SBATCH --job-name="theo_spec"
#SBATCH --output=theo_spectra.out
#SBATCH --partition=accel
#SBATCH --time=30-00:00:00

module purge
module load system/R−3.6.1

~/R-3.6.1/bin/Rscript sim1.R
