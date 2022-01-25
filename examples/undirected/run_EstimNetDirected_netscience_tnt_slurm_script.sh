#!/bin/bash

#SBATCH --job-name="tnt_EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-0:10:00
#SBATCH --output=EstimNetDirected_tnt_netscience-%j.out
#SBATCH --error=EstimNetDirected_tnt_netscience-%j.err

echo -n "started at: "; date

ROOT=../..

module purge
module load openmpi

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_netscience_tnt.txt

module load r

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_tnt_netscience dzA_tnt_netscience
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_tnt_netscience dzA_tnt_netscience

echo -n "ended at: "; date

