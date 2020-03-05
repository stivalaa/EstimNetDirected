#!/bin/bash

#SBATCH --job-name="SimualateERGM_GoF"
#SBATCH --ntasks=1
#SBATCH --time=0-4:00:00

echo -n "started at: "; date

ROOT=..

module load R/3.2.5

time ${ROOT}/src/SimulateERGM  sim_config_example_estimated_contattr.txt

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_estimated_sim_n2000_cont_sample.txt  obs_stats_n2000_cont_sample_0.txt

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R  sample_statistics_n2000_directed_cont_sim7920000000.txt simulation_estimated_sim_n2000_cont_sample

echo -n "ended at: "; date

