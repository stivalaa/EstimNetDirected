#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi_notables"
#SBATCH --ntasks=8
#SBATCH --time=0-2:00:00
#SBATCH --output=EstimNetDirected_notables_contattr_example-%j.out
#SBATCH --error=EstimNetDirected_notables_contattr_example-%j.err

echo -n "started at: "; date

ROOT=..

module load openmpi

time mpirun ${ROOT}/src/EstimNetDirected_mpi config_example_contattr.txt

echo -n "ended at: "; date

