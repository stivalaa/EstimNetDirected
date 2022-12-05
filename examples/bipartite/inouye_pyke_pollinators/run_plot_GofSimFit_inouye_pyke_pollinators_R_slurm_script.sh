#!/bin/bash

#SBATCH --job-name="plotGofSimFit-inouye_pyke_pollinators"
#SBATCH --ntasks=1
#SBATCH --time=0-00:10:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --partition=slim
#SBATCH --output=PlotGofSimFit-inouye_pyke_pollinators-%j.out
#SBATCH --error=PlotGofSimFit-inouye_pyke_pollinators-%j.err

echo -n "started at: "; date

ROOT=${HOME}/EstimNetDirected

module load r

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R -y -t  inouye_pyke_pollinators_bipartite.net sim_gof_inouye_pyke_pollinators

echo -n "ended at: "; date

