#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi_hash"
#SBATCH --ntasks=8
#SBATCH --time=0-1:00:00

echo -n "started at: "; date

ROOT=..

module purge
module load openmpi

time mpirun ${ROOT}/src/EstimNetDirected_mpi_hashtables config_polblogs.txt

module load r

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_polblogs dzA_polblogs
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_polblogs dzA_polblogs

echo -n "ended at: "; date

