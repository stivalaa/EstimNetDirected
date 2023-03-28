#!/bin/bash

#SBATCH --job-name="plotEstimNetDirectedSimFit-robertson_pollinators"
#SBATCH --ntasks=1
#SBATCH --time=0-1:00:00
#SBATCH --mem-per-cpu=6GB
#SBATCH --partition=slim
#SBATCH --output=PlotEstimNetDirectedSimFit-robertson_pollinators-%j.out
#SBATCH --error=PlotEstimNetDirectedSimFit-robertson_pollinators-%j.err

echo -n "started at: "; date

ROOT=../../..

module load r

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R ../../../Test/TestChangeStatsBipartite/robertson_pollinators_bipartite.net sim_ifd_robertson_pollinators

echo -n "ended at: "; date

