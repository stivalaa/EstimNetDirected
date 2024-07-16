#!/bin/bash

#SBATCH --job-name="SimulateERGM-inouye_pyke_pollinators_altkcycles"
#SBATCH --ntasks=1
#SBATCH --time=0-0:10:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=SimulateERGM-inouye_pyke_pollinators_altkcycles-%j.out
#SBATCH --error=SimulateERGM-inouye_pyke_pollinators_altkcycles-%j.err

echo -n "started at: "; date
uname -a

#ROOT=${HOME}/EstimNetDirected
ROOT=../../../

module purge
module load gcc/11.3.0    # need this for 'toolchain' hierarchical modules
module load openmpi/4.1.4 # module version numbers are required on OzStar

gof_sim_config=config_gof_inouye_pyke_pollinators_altkcycles.txt

if [ -f ${gof_sim_config} ]; then
  mv ${gof_sim_config} ${gof_sim_config}.OLD
fi

rm sim_gof_inouye_pyke_pollinators_altkcycles_*.net

${ROOT}/scripts/estimnetdirectedEstimation2simulationConfig.sh config_inouye_pyke_pollinators_ifd_altkcycles.txt estimation_ifd_inouye_pyke_pollinators_altkcycles.txt stats_gof_inouye_pyke_pollinators_altkcycles.txt sim_gof_inouye_pyke_pollinators_altkcycles  | sed 's/interval = 100000 /interval = 10000 /' | sed 's/burnin = 1000000 /burnin = 100000 /' | sed 's/TNT/IFD/g'  > ${gof_sim_config}
echo 'numArcs = 281' >> ${gof_sim_config}


time ${ROOT}/src/SimulateERGM ${gof_sim_config}

module load r/4.2.1

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_gof_inouye_pyke_pollinators_altkcycles.txt  obs_stats_ifd_inouye_pyke_pollinators_altkcycles_0.txt

echo -n "ended at: "; date

