#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi_hashtables"
#SBATCH --ntasks=8
#SBATCH --time=0-2:00:00
#SBATCH --output=EstimNetDirected_hashtables_contattr_example-%j.out
#SBATCH --error=EstimNetDirected_hashtables_contattr_example-%j.err

echo -n "started at: "; date

ROOT=..

module load openmpi

time mpirun ${ROOT}/src/EstimNetDirected_mpi_hashtables config_example_contattr.txt

echo -n "ended at: "; date

