#!/bin/bash

#SBATCH --job-name="SimulateERGM"
#SBATCH --ntasks=1
#SBATCH --time=0-2:30:00

echo -n "started at: "; date

ROOT=..

time ${ROOT}/src/SimulateERGM  sim_config_example.txt

module load r

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_sim_n2000_sample.txt

echo -n "ended at: "; date

