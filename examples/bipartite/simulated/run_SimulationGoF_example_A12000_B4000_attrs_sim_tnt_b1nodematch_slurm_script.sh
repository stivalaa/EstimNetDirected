#!/bin/bash

#SBATCH --job-name="SimulationGoF_mpi_tnt_b1nodematch_bpnet_A12000_B4000_sim"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=0-1:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --partition=slim
#SBATCH --output=SimulationGoF_tnt_b1nodematch_bpnet_A12000_B4000_sim-%j.out
#SBATCH --error=SimulationGoF_tnt_b1nodematch_bpnet_A12000_B4000_sim-%j.err

echo -n "started at: "; date

ROOT=../../..

gof_sim_config=config_gof_bipartite_A12000_B4000_attrs_sim_tnt_b1nodematch.txt

if [ -f ${gof_sim_config} ]; then
  mv ${gof_sim_config} ${gof_sim_config}.OLD
fi

rm sim_gof_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch_sim_*.net

module purge

${ROOT}/scripts/estimnetdirectedEstimation2simulationConfig.sh  config_bipartite_A12000_B4000_attrs_sim_tnt_b1nodematch.txt estimation_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch.txt stats_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch_sim.txt sim_gof_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch_sim  > ${gof_sim_config}

time ${ROOT}/src/SimulateERGM ${gof_sim_config}

module load r 

# for R library(parallel) mclapply()
export MC_CORES=${SLURM_CPUS_ON_NODE}
echo MC_CORES = $MC_CORES

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R  stats_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch_sim.txt obs_stats_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch_0.txt

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R -y 16 -s  bpnet_A12000_B4000_attrs_sim830000000.net sim_gof_tnt_bpnet_A12000_B4000_attrs_sim_b1nodematch_sim


echo -n "ended at: "; date

