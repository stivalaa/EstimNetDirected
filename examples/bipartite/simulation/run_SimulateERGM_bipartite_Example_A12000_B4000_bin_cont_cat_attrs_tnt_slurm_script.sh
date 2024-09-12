#!/bin/bash

#SBATCH --job-name="SimulateERGM"
#SBATCH --ntasks=1
#SBATCH --time=0-0:30:00

echo -n "started at: "; date
uname -a

ROOT=../../..

module load r

time ${ROOT}/src/SimulateERGM  config_sim_bipartite_example_A12000_B4000_bin_cont_cat_attrs_tnt.txt

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_sim_bipartite_A12000_B4000_bin_cont_cat_attrs_sampler_tnt.txt

echo -n "ended at: "; date

