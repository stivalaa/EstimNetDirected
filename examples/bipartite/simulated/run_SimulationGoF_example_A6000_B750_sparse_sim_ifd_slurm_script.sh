#!/bin/bash

#SBATCH --job-name="SimulationGoF_mpi_ifd_bpnet_A6000_B750_sim"
#SBATCH --ntasks=1
#SBATCH --time=0-0:20:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=SimulationGoF_ifd_bpnet_A6000_B750_sim-%j.out
#SBATCH --error=SimulationGoF_ifd_bpnet_A6000_B750_sim-%j.err

echo -n "started at: "; date

ROOT=../../..

gof_sim_config=config_gof_ifd_bpnet_A6000_B750_sim.txt

if [ -f ${gof_sim_config} ]; then
  mv ${gof_sim_config} ${gof_sim_config}.OLD
fi

rm sim_gof_ifd_bpnet_A6000_B750_sim_*.net

module purge

${ROOT}/scripts/estimnetdirectedEstimation2simulationConfig.sh config_bipartite_A6000_B750_sparse_sim_ifd.txt  estimation_ifd_bpnet_A6000_B750_sim.txt  stats_ifd_bpnet_A6000_B750_sim.txt sim_gof_ifd_bpnet_A6000_B750_sim > ${gof_sim_config}

time ${ROOT}/src/SimulateERGM ${gof_sim_config}

module load r 

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_ifd_bpnet_A6000_B750_sim.txt obs_stats_ifd_bpnet_A6000_B750_sim_0.txt 

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R bpnet_A6000_B750_sparse_sim100000000.net sim_gof_ifd_bpnet_A6000_B750_sim


echo -n "ended at: "; date

