#!/bin/bash

#SBATCH --job-name="plotEstimNetDirectedSimFit-inouye_pyke_pollinators"
#SBATCH --ntasks=1
#SBATCH --time=0-1:00:00
#SBATCH --mem-per-cpu=6GB
#SBATCH --partition=slim
#SBATCH --output=PlotEstimNetDirectedSimFit-inouye_pyke_pollinators-%j.out
#SBATCH --error=PlotEstimNetDirectedSimFit-inouye_pyke_pollinators-%j.err

echo -n "started at: "; date

#ROOT=${HOME}/EstimNetDirected
ROOT=../../../

module load r

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R -y 16 inouye_pyke_pollinators_bipartite.net sim_ifd_inouye_pyke_pollinators

echo -n "ended at: "; date

