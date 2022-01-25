#!/bin/bash

#SBATCH --job-name="basic_EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-0:10:00
#SBATCH --output=EstimNetDirected_basic_netscience-%j.out
#SBATCH --error=EstimNetDirected_basic_netscience-%j.err

echo -n "started at: "; date

ROOT=../..

module purge
module load openmpi

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_netscience_basic.txt

module load r

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_basic_netscience dzA_basic_netscience
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_basic_netscience dzA_basic_netscience

echo -n "ended at: "; date

