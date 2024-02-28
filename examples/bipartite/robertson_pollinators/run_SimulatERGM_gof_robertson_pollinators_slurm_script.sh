#!/bin/bash

#SBATCH --job-name="SimulateERGM-robertson_pollinators"
#SBATCH --ntasks=1
#SBATCH --time=0-24:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=SimulateERGM-robertson_pollinators-%j.out
#SBATCH --error=SimulateERGM-robertson_pollinators-%j.err

echo -n "started at: "; date
uname -a

ROOT=../../..

module purge
module load gcc/11.3.0    # need this for 'toolchain' hierarchical modules
module load openmpi/4.1.4 # module version numbers are required on Oz

gof_sim_config=config_gof_robertson_pollinators.txt

if [ -f ${gof_sim_config} ]; then
  mv ${gof_sim_config} ${gof_sim_config}.OLD
fi

rm sim_gof_robertson_pollinators_*.net

${ROOT}/scripts/estimnetdirectedEstimation2simulationConfig.sh config_robertson_pollinators_ifd.txt estimation_ifd_robertson_pollinators.txt stats_gof_robertson_pollinators.txt sim_gof_robertson_pollinators  | sed 's/interval = 100000 /interval = 10000 /' | sed 's/burnin = 1000000 /burnin = 100000 /' | sed 's/TNT/IFD/g'  > ${gof_sim_config}
echo 'numArcs = 15255' >> ${gof_sim_config}


time ${ROOT}/src/SimulateERGM ${gof_sim_config}

module load r/4.2.1

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_gof_robertson_pollinators.txt  obs_stats_ifd_robertson_pollinators_0.txt

echo -n "ended at: "; date

