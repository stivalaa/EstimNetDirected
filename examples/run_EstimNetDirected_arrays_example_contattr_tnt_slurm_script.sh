#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi_arrays"
#SBATCH --ntasks=8
#SBATCH --time=0-2:00:00
#SBATCH --output=EstimNetDirected_arrays_contattr_example-%j.out
#SBATCH --error=EstimNetDirected_arrays_contattr_example-%j.err

echo -n "started at: "; date

ROOT=..

module load openmpi

time mpirun ${ROOT}/src/EstimNetDirected_mpi_arrays config_example_contattr.txt

echo -n "ended at: "; date

