#!/bin/bash

#SBATCH --job-name="plotEstimNetDirectedSimFit-robertson_pollinators"
#SBATCH --ntasks=1
#SBATCH --time=0-1:00:00
#SBATCH --mem-per-cpu=6GB
#SBATCH --output=PlotEstimNetDirectedSimFit-robertson_pollinators-%j.out
#SBATCH --error=PlotEstimNetDirectedSimFit-robertson_pollinators-%j.err

echo -n "started at: "; date

ROOT=../../..

command -v module > /dev/null 2>&1 && module load gcc/11.3.0 # needed by r/4.2.1
command -v module > /dev/null 2>&1 && module load openmpi/4.1.4 # needed by r/4.2.1
command -v module > /dev/null 2>&1 && module load r/4.2.1

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R ../../../Test/TestChangeStatsBipartite/robertson_pollinators_bipartite.net sim_ifd_robertson_pollinators

echo -n "ended at: "; date

