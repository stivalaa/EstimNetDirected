#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi_basic_example_n1000_binattr"
#SBATCH --ntasks=32
#SBATCH --time=0-0:30:00
#SBATCH --output=EstimNetDirected_basic_example_n1000_binattr-%j.out
#SBATCH --error=EstimNetDirected_basic_example_n1000_binattr-%j.err

echo -n "started at: "; date

ROOT=../..

module load openmpi

#cannot do this here on upgraded system as loding R module causes many things (and specifically anything using MPI) to fail: module load r

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_example_n1000_binattar_basic.txt

module load r # have to do this AFTER running MPI progams on 'upgraded' cluster (also not 'r' not 'R')

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_basic_example_n1000_binattr dzA_basic_example_n1000_binattr | tee estimation_basic_example_n1000_binattr.txt
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_basic_example_n1000_binattr dzA_basic_example_n1000_binattr
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R sample_statistics_n1000_binattr_50_50_sim272500000.txt  sim_basic_example_n1000_binattr

echo -n "ended at: "; date

