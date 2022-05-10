#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi_basic_bpnet_A2000_B250_sim"
#SBATCH --ntasks=32
#SBATCH --time=0-2:00:00
#SBATCH --output=EstimNetDirected_basic_bpnet_A2000_B250_sim-%j.out
#SBATCH --error=EstimNetDirected_basic_bpnet_A2000_B250_sim-%j.err

echo -n "started at: "; date

ROOT=../../..

module purge
module load openmpi

#cannot do this here on upgraded system as loding R module causes many things (and specifically anything using MPI) to fail: module load r

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_bipartite_A750_B250_sim_basic.txt

module load r # have to do this AFTER running MPI progams on 'upgraded' cluster (also not 'r' not 'R')

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_basic_bpnet_A2000_B250_sim dzA_basic_bpnet_A2000_B250_sim | tee estimation_basic_bpnet_A2000_B250_sim.txt
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_basic_bpnet_A2000_B250_sim dzA_basic_bpnet_A2000_B250_sim
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R bpnet_A2000_B250_sim5450000000.net  sim_basic_bpnet_A2000_B250_sim

echo -n "ended at: "; date

