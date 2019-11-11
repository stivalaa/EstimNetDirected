#!/bin/bash

#SBATCH --job-name="tnt_EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-0:30:00

echo -n "started at: "; date

ROOT=..

module load openmpi
module load R

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_example_tnt.txt

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_tnt_sim_n1000_sample dzA_tnt_sim_n1000_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_tnt_sim_n1000_sample dzA_tnt_sim_n1000_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R ../pythonDemo/sample_statistics_n1000_directed_binattr_sim620000000.txt sim_tnt_sim_n1000_sample

echo -n "ended at: "; date

