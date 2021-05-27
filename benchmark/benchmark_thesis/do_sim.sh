#!/bin/bash
#SBATCH --output=benchmark.out
#SBATCH --job-name=benchmark
#SBATCH --partition=accel
#SBATCH --time=30-00:00:00
module purge
module load system/R−3.6.1

~/R-3.6.1/bin/Rscript benchmark.R