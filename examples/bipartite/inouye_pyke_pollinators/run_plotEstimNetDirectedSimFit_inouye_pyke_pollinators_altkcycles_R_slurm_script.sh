#!/bin/bash

#SBATCH --job-name="plotEstimNetDirectedSimFit-inouye_pyke_pollinators_altkcycles"
#SBATCH --ntasks=1
#SBATCH --time=0-1:00:00
#SBATCH --mem-per-cpu=6GB
#SBATCH --output=PlotEstimNetDirectedSimFit-inouye_pyke_pollinators_altkcycles-%j.out
#SBATCH --error=PlotEstimNetDirectedSimFit-inouye_pyke_pollinators_altkcycles-%j.err

echo -n "started at: "; date

#ROOT=${HOME}/EstimNetDirected
ROOT=../../../

command -v module > /dev/null 2>&1 && module load gcc/11.3.0 # needed by r/4.2.1
command -v module > /dev/null 2>&1 && module load openmpi/4.1.4 # needed by r/4.2.1
command -v module > /dev/null 2>&1 && module load r/4.2.1

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R -y 6 -t inouye_pyke_pollinators_altkcycles_bipartite.net sim_ifd_inouye_pyke_pollinators_altkcycles

## too slow to count cycles (up to 16):
#time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R -y 16 inouye_pyke_pollinators_altkcycles_bipartite.net sim_ifd_inouye_pyke_pollinators_altkcycles

echo -n "ended at: "; date

