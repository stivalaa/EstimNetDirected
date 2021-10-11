#!/bin/bash

#SBATCH --job-name="ifd_EstimNetDirected_mpi"
#SBATCH --ntasks=32
#SBATCH --time=0-0:30:00

echo -n "started at: "; date

ROOT=..

module purge
module load openmpi

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_example_ifd_allowloops.txt

module load r

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_ifd_allowloops_sim_n1000_sample dzA_ifd_allowloops_sim_n1000_sample | tee estimation_ifd_allowloops_sim_n1000_sample.out
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_ifd_allowloops_sim_n1000_sample dzA_ifd_allowloops_sim_n1000_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R ../pythonDemo/sample_statistics_n1000_directed_binattr_sim620000000.txt sim_ifd_allowloops_sim_n1000_sample

echo -n "ended at: "; date

