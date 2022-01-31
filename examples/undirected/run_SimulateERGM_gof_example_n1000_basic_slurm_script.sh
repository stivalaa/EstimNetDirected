#!/bin/bash

#SBATCH --job-name="SimulateERGM_gof_example_n1000_binattr"
#SBATCH --ntasks=1
#SBATCH --time=0-0:30:00
#SBATCH --output=SimulateERGM_gof_basic_example_n1000_binattr-%j.out
#SBATCH --error=SimulateERGM_gof_basic_example_n1000_binattr-%j.err

echo -n "started at: "; date

ROOT=../..

SIM_FILE_PREFIX=sim_gof_example_n1000_binattr_basic
STATS_FILE=stats_gof_example_n1000_binattr_basic.txt
GOF_CONFIG_FILE=config_gof_example_n1000_binattr_basic.txt

module load r

${ROOT}/scripts/estimnetdirectedEstimation2simulationConfig.sh  config_example_n1000_binattar_basic.txt  estimation_basic_example_n1000_binattr.txt ${STATS_FILE} ${SIM_FILE_PREFIX} > ${GOF_CONFIG_FILE}

rm ${SIM_FILE_PREFIX}_*.net

time ${ROOT}/src/SimulateERGM  ${GOF_CONFIG_FILE}

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R  ${STATS_FILE}  obs_stats_basic_example_n1000_binattr_0.txt

time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R  sample_statistics_n1000_binattr_50_50_sim272500000.txt ${SIM_FILE_PREFIX}


echo -n "ended at: "; date

