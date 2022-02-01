#!/bin/bash

#SBATCH --job-name="SimulateERGM_gof_netscience"
#SBATCH --ntasks=1
#SBATCH --time=0-0:30:00
#SBATCH --output=SimulateERGM_gof_tnt_netscience-%j.out
#SBATCH --error=SimulateERGM_gof_tnt_netscience-%j.err

echo -n "started at: "; date

ROOT=../..

SIM_FILE_PREFIX=sim_gof_netscience_tnt
STATS_FILE=stats_gof_netscience_tnt.txt
GOF_CONFIG_FILE=config_gof_netscience_tnt.txt

module load r

${ROOT}/scripts/estimnetdirectedEstimation2simulationConfig.sh  config_netscience_tnt.txt   estimation_tnt_netscience.txt ${STATS_FILE} ${SIM_FILE_PREFIX} > ${GOF_CONFIG_FILE}

rm ${SIM_FILE_PREFIX}_*.net

time ${ROOT}/src/SimulateERGM  ${GOF_CONFIG_FILE}

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R  ${STATS_FILE}  obs_stats_tnt_netscience_0.txt

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R  ../../Test/TestChangeStatsUndirected/netscience_edgelist.txt ${SIM_FILE_PREFIX}


echo -n "ended at: "; date

