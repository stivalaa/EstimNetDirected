#!/bin/bash

#SBATCH --job-name="ifd_b_EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-1:30:00

echo -n "started at: "; date

ROOT=..

module load openmpi
module load R

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_polblogs_ifd_borisenko.txt

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_ifd_borisenko_polblogs dzA_ifd_borisenko_polblogs
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_ifd_borisenko_polblogs dzA_ifd_borisenko_polblogs

echo -n "ended at: "; date

