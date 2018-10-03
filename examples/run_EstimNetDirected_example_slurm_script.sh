#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-0:30:00

echo -n "started at: "; date

ROOT=..

module load openmpi
module load R

time srun ${ROOT}/src/EstimNetDirected_mpi config_example.txt

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_sim_n1000_sample dzA_sim_n1000_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_sim_n1000_sample dzA_sim_n1000_sample

echo -n "ended at: "; date

