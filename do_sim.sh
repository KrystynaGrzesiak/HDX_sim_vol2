#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --job-name="A long job"
#SBATCH --mem=5GB
#SBATCH --output=long-job.out

module purge

module load system/R−3.6.1

~/R-3.6.1/bin/Rscript theoretical_spectra_sim.R
