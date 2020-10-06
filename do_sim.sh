#!/bin/bash
#SBATCH --job-name="theo_spec"
#SBATCH --output=theo_spectra.out

module purge

module load system/R−3.6.1

~/R-3.6.1/bin/Rscript theoretical_spectra_sim.R
