#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-6:00:00

echo -n "started at: "; date

ROOT=..

module load R/3.2.1-vlsci_intel-2015.08.25

time srun ${ROOT}/src/EstimNetDirected_mpi config_polblogs.txt

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_polblogs dzA_polblogs
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_polblogs dzA_polblogs

echo -n "ended at: "; date

