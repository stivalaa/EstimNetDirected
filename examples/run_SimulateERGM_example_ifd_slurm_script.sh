#!/bin/bash

#SBATCH --job-name="SimualateERGM_IFD"
#SBATCH --ntasks=1
#SBATCH --time=0-0:30:00

echo -n "started at: "; date

ROOT=..

module load R

time ${ROOT}/src/SimulateERGM  sim_config_example_ifd.txt

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_sim_n2000_sample_ifd.txt

echo -n "ended at: "; date

