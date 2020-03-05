#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-0:30:00

echo -n "started at: "; date

ROOT=..

module load openmpi
module load R/3.2.5

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_example_contattr.txt

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_sim_n2000_cont_sample dzA_sim_n2000_cont_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_sim_n2000_cont_sample dzA_sim_n2000_cont_sample

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R sample_statistics_n2000_directed_cont_sim7920000000.txt sim_sim_n2000_cont_sample

echo -n "ended at: "; date

