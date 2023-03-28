#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi-ifd_robertson_pollinators"
#SBATCH --ntasks=20
#SBATCH --time=0-00:10:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --partition=slim
#SBATCH --output=EstimNetDirected-ifd_robertson_pollinators-%j.out
#SBATCH --error=EstimNetDirected-ifd_robertson_pollinators-%j.err

echo -n "started at: "; date

ROOT=../../../

module purge
module load openmpi


cp -p ${ROOT}/src/EstimNetDirected_mpi .

time mpirun ./EstimNetDirected_mpi config_robertson_pollinators_ifd.txt

module load r

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_ifd_robertson_pollinators dzA_ifd_robertson_pollinators |tee estimation_ifd_robertson_pollinators.txt
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_ifd_robertson_pollinators dzA_ifd_robertson_pollinators

echo -n "ended at: "; date

