#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi"
#SBATCH --ntasks=32
#SBATCH --time=0-0:10:00
#SBATCH --mem-per-cpu=1G

echo -n "started at: "; date

ROOT=..

module load gcc/9.2.0     # need this for 'toolchain' hierarchical modules
module load openmpi/4.0.2 # module version numbers are required on OzStar

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_example.txt

module load r/4.1.1 # module version numbers are required on OzStar

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_sim_n1000_sample dzA_sim_n1000_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_sim_n1000_sample dzA_sim_n1000_sample
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R ../pythonDemo/sample_statistics_n1000_directed_binattr_sim620000000.txt sim_sim_n1000_sample

echo -n "ended at: "; date

