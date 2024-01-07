#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch"
#SBATCH --ntasks=64
#SBATCH --time=0-1:00:00
#SBATCH --mem-per-cpu=512MB
#SBATCH --output=EstimNetDirected_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch-%j.out
#SBATCH --error=EstimNetDirected_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch-%j.err

echo -n "started at: "; date

ROOT=../../..

module purge
module load openmpi

#cannot do this here on upgraded system as loding R module causes many things (and specifically anything using MPI) to fail: module load r

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_bipartite_A12000_B4000_attrs_sim_tnt_b1nodematch.txt

module load r # have to do this AFTER running MPI progams on 'upgraded' cluster (also not 'r' not 'R')

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch dzA_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch | tee estimation_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch.txt
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch dzA_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R bpnet_A12000_B4000_attrs_sim830000000.net sim_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch

echo -n "ended at: "; date

