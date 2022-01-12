#!/bin/bash

#SBATCH --job-name="ifd_b_EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-6:00:00

echo -n "started at: "; date

ROOT=..

module purge
module load openmpi

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_polblogs_ifd_borisenko.txt

module load r

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_ifd_borisenko_polblogs dzA_ifd_borisenko_polblogs
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_ifd_borisenko_polblogs dzA_ifd_borisenko_polblogs

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R ../pythonDemo/polblogs/polblogs_arclist.txt sim_ifd_borisenko_polblogs

echo -n "ended at: "; date

