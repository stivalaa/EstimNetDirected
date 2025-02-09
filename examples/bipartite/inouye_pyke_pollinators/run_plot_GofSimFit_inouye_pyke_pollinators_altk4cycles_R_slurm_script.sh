#!/bin/bash

#SBATCH --job-name="plotGofSimFit-inouye_pyke_pollinators_altk4cycles"
#SBATCH --ntasks=1
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1GB
#SBATCH --output=PlotGofSimFit-inouye_pyke_pollinators_altk4cycles-%j.out
#SBATCH --error=PlotGofSimFit-inouye_pyke_pollinators_altk4cycles-%j.err

echo -n "started at: "; date

#ROOT=${HOME}/EstimNetDirected
ROOT=../../../

module load gcc/11.3.0 # needed by r/4.2.1
module load openmpi/4.1.4 # needed by r/4.2.1
module load r/4.2.1

# for R library(parallel) mclapply()
export MC_CORES=${SLURM_CPUS_ON_NODE}
echo MC_CORES = $MC_CORES

## far too slow with -y 16 (obs not done afer 25 mins)
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R -y 12 -t  -s inouye_pyke_pollinators_bipartite.net sim_gof_inouye_pyke_pollinators_altk4cycles

echo -n "ended at: "; date

