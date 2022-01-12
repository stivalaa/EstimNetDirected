#!/bin/bash

#SBATCH --job-name="SimualateERGM"
#SBATCH --ntasks=1
#SBATCH --time=0-4:00:00

echo -n "started at: "; date

ROOT=..

module load r

time ${ROOT}/src/SimulateERGM  sim_config_example_binattr.txt

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_sim_n1000_binattr_sample.txt

echo -n "ended at: "; date

