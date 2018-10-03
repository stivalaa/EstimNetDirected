#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-0:30:00

echo -n "started at: "; date

ROOT=..

module load openmpi
module load R

time srun ${ROOT}/src/EstimNetDirected_mpi config_polblogs.txt

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_polblogs dzA_polblogs
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_polblogs dzA_polblogs

echo -n "ended at: "; date

