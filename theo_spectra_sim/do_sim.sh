#!/bin/bash
#SBATCH --job-name="theo_spec"
#SBATCH --output=theo_spectra.out
#SBATCH --partition=accel
#SBATCH --time=30-00:00:00

module purge
module load ~/R/
source ~/setup.sh %
Rscript  ~/HDX_sim_vol2/theo_spectra_sim/theoretical_spectra_sim.R




