#!/bin/bash

#SBATCH --job-name="SimualateERGM"
#SBATCH --ntasks=1
#SBATCH --time=0-2:30:00

echo -n "started at: "; date

ROOT=..

module load R

time ${ROOT}/src/SimulateERGM  sim_config_estimated_polblogs.txt

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_sim_estimated_polblogs.txt  obs_stats_ifd_boriseko_polblogs_0.txt

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R  ../pythonDemo/polblogs/polblogs_arclist.txt simulation_estimated_polblogs


echo -n "ended at: "; date

