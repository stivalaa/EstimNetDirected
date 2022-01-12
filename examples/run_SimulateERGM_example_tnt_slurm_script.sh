#!/bin/bash

#SBATCH --job-name="SimulateERGM"
#SBATCH --ntasks=1
#SBATCH --time=0-2:30:00

echo -n "started at: "; date

ROOT=..

module load r

time ${ROOT}/src/SimulateERGM  sim_config_example_tnt.txt

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_sim_n2000_sample_tnt.txt

echo -n "ended at: "; date

