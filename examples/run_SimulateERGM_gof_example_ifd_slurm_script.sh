#!/bin/bash

#SBATCH --job-name="SimualateERGM_IFD_GOF"
#SBATCH --ntasks=1
#SBATCH --time=0-3:30:00

echo -n "started at: "; date

ROOT=..


config_tmpfile=sim_config_estimated_ifd_sim_n1000_sample.txt
statsfile=stats_estimated_ifd_sim_n1000_sample.txt
simfileprefix=gof_estimated_ifd_sim_n1000_sample


rm ${simfileprefix}_*.net

# generate SimulationERGM config file with parameter values from estimation
${ROOT}/scripts/estimnetdirectedEstimation2simulationConfig.sh  config_example_ifd.txt estimation_ifd_sim_n1000_sample.out ${statsfile} ${simfileprefix} > ${config_tmpfile}

# simulate networks from the estimated parameters
time ${ROOT}/src/SimulateERGM ${config_tmpfile}

module load r

# plot simulation diagnostics with observed sufficient statistics 
time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R ${statsfile} obs_stats_ifd_n1000_sample_0.txt

# plot goodness-of-fit observed and simulated network statistics (not in model)
time Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R ../pythonDemo/sample_statistics_n1000_directed_binattr_sim620000000.txt ${simfileprefix}


echo -n "ended at: "; date

