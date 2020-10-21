#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi"
#SBATCH --ntasks=32
#SBATCH --time=0-0:30:00

echo -n "started at: "; date

ROOT=..

module load openmpi

#cannot do this here on upgraded system as loding R module causes many things (and specifically anything using MPI) to fail: module load r

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_example.txt

module load r # have to do this AFTER running MPI progams on 'upgraded' cluster (also not 'r' not 'R')

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_sim_n1000_sample dzA_sim_n1000_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_sim_n1000_sample dzA_sim_n1000_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R ../pythonDemo/sample_statistics_n1000_directed_binattr_sim620000000.txt sim_sim_n1000_sample

echo -n "ended at: "; date

