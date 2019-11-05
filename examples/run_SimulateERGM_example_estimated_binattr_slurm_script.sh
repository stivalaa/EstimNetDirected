#!/bin/bash

#SBATCH --job-name="SimualateERGM_GoF"
#SBATCH --ntasks=1
#SBATCH --time=0-4:00:00

echo -n "started at: "; date

ROOT=..

module load R

time ${ROOT}/src/SimulateERGM  sim_config_example_estimated_binattr.txt

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_estimated_sim_n1000_binattr_sample.txt  obs_stats_n1000_sample_0.txt

echo -n "ended at: "; date

