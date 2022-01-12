#!/bin/bash

#SBATCH --job-name="SimualateERGM_GoF"
#SBATCH --ntasks=1
#SBATCH --time=0-4:00:00

echo -n "started at: "; date

ROOT=..

module load r

time ${ROOT}/src/SimulateERGM  sim_config_example_estimated_binattr.txt

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_estimated_sim_n1000_binattr_sample.txt  obs_stats_n1000_sample_0.txt

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R  ../pythonDemo/sample_statistics_n1000_directed_binattr_sim620000000.txt simulation_estimated_sim_n1000_binattr_sample

echo -n "ended at: "; date

