#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi"
#SBATCH --ntasks=32
#SBATCH --time=0-0:10:00
#SBATCH --mem-per-cpu=1G

echo -n "started at: "; date

ROOT=..

module load foss/2022b
module load python/3.10.8
module load numpy/1.24.2-scipy-bundle-2023.02

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_example.txt

module load gcc/11.3.0    # needed by r/4.2.1
module load openmpi/4.1.4 # needed by r/4.2.1
module load r/4.2.1

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_sim_n1000_sample dzA_sim_n1000_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_sim_n1000_sample dzA_sim_n1000_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R ../pythonDemo/sample_statistics_n1000_directed_binattr_sim620000000.txt sim_sim_n1000_sample

echo -n "ended at: "; date

