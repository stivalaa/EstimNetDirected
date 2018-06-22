#!/bin/bash

#SBATCH --job-name="ifd_EstimNetDirected_mpi"
#SBATCH --ntasks=8
#SBATCH --time=0-1:30:00

echo -n "started at: "; date

ROOT=..

module load R/3.2.1-vlsci_intel-2015.08.25

time srun ${ROOT}/src/EstimNetDirected_mpi config_polblogs_ifd.txt

time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_ifd_polblogs dzA_ifd_polblogs
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_ifd_polblogs dzA_ifd_polblogs

echo -n "ended at: "; date

