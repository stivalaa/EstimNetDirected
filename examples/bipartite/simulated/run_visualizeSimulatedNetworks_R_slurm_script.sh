#!/bin/bash

#SBATCH --job-name="visualizeSimulatedNetworks.R"
#SBATCH --ntasks=1
#SBATCH --time=0-8:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --partition=slim
#SBATCH --output=visualizeSimulatedNetworks-%j.out
#SBATCH --error=visualizeSimulatedNetworks-%j.err

echo -n "started at: "; date

module purge
module load r

time Rscript visualizeSimulatedGraphs.R

echo -n "ended at: "; date

