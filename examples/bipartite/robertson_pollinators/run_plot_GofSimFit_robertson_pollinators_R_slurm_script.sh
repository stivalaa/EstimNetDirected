#!/bin/bash

#SBATCH --job-name="plotGofSimFit-robertson_pollinators"
#SBATCH --ntasks=1
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1GB
#SBATCH --partition=slim
#SBATCH --output=PlotGofSimFit-robertson_pollinators-%j.out
#SBATCH --error=PlotGofSimFit-robertson_pollinators-%j.err

echo -n "started at: "; date

ROOT=../../..

module load r

# for R library(parallel) mclapply()
export MC_CORES=${SLURM_CPUS_ON_NODE}
echo MC_CORES = $MC_CORES

# even -y 8 does not complete in an hour (started at 12 and worked down)
# and -y 6 over half an hour
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R  -s ../../../Test/TestChangeStatsBipartite/robertson_pollinators_bipartite.net sim_gof_robertson_pollinators

echo -n "ended at: "; date
