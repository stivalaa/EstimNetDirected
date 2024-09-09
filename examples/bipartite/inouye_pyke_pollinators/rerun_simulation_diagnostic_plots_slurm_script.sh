#!/bin/bash

#SBATCH --job-name="rerun_diagnostic_plots_inouye_pyke_pollinators"
#SBATCH --ntasks=1
#SBATCH --time=0-0:10:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=rerun_diagnostic_plots_inouye_pyke_pollinators-%j.out
#SBATCH --error=rerun_diagnostic_plots_inouye_pyke_pollinators-%j.err

echo -n "started at: "; date
uname -a

#ROOT=${HOME}/EstimNetDirected
ROOT=../../../

module purge
module load gcc/11.3.0    # need this for 'toolchain' hierarchical modules
module load openmpi/4.1.4 # module version numbers are required on OzStar
module load r/4.2.1

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_gof_inouye_pyke_pollinators.txt  obs_stats_ifd_inouye_pyke_pollinators_0.txt

time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R stats_gof_inouye_pyke_pollinators_altkcycles.txt  obs_stats_ifd_inouye_pyke_pollinators_altkcycles_0.txt

echo -n "ended at: "; date

