#!/bin/bash

#SBATCH --job-name="R_contattr_Example"
#SBATCH --ntasks=1
#SBATCH --time=0-0:10:00
#SBATCH --output=R_contattr_example_EstimNetDirected_estim_simfit-%j.out
#SBATCH --error=R_contattr_example_EstimNetDirected_estim_simfit-%j.err

echo -n "started at: "; date

ROOT=..

module load R/3.2.5


time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_sim_n2000_cont_sample dzA_sim_n2000_cont_sample

time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_sim_n2000_cont_sample dzA_sim_n2000_cont_sample

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R sample_statistics_n2000_directed_cont_sim7920000000.txt sim_sim_n2000_cont_sample

echo -n "ended at: "; date

