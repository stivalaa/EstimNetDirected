#!/bin/bash

#SBATCH --job-name="tnt_EstimNetDirected_mpi_hippie"
#SBATCH --ntasks=32
#SBATCH --time=0-2:00:00
#SBATCH --output=EstimNetDirected_tnt_hippie-%j.out
#SBATCH --error=EstimNetDirected_tnt_hippie-%j.err

echo -n "started at: "; date

ROOT=../..

module purge
module load openmpi

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_hippie_tnt.txt

module load r

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_tnt_hippie dzA_tnt_hippie | tee estimation_tnt_hippie.txt
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_tnt_hippie dzA_tnt_hippie
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R ../../examples/undirected/hippie_ppi_high_edgelist.txt sim_tnt_hippie

echo -n "ended at: "; date

