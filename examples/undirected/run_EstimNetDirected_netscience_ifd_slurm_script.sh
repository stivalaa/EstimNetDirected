#!/bin/bash

#SBATCH --job-name="ifd_EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-0:10:00
#SBATCH --output=EstimNetDirected_ifd_netscience-%j.out
#SBATCH --error=EstimNetDirected_ifd_netscience-%j.err

echo -n "started at: "; date

ROOT=../..

module purge
module load openmpi

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_netscience_ifd.txt

module load r

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_ifd_netscience dzA_ifd_netscience
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_ifd_netscience dzA_ifd_netscience

echo -n "ended at: "; date

