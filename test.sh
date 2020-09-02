#!/bin/bash
#SBATCH −J lauchRscript
#SBATCH −ooutput . out

module purge

module load system/R−3.5.1

Rscript test.R