#!/bin/bash

#SBATCH --job-name="ifd_EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-0:30:00

echo -n "started at: "; date

ROOT=..

module load openmpi
module load R

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_example_ifd.txt

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_ifd_sim_n1000_sample dzA_ifd_sim_n1000_sample | tee estimation_ifd_sim_n1000_sample.out
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_ifd_sim_n1000_sample dzA_ifd_sim_n1000_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R ../pythonDemo/sample_statistics_n1000_directed_binattr_sim620000000.txt sim_ifd_sim_n1000_sample

echo -n "ended at: "; date

