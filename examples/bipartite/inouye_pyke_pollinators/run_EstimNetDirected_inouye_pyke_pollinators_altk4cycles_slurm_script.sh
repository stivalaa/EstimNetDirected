#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi-ifd_inouye_pyke_pollinators_altk4cycles"
#SBATCH --ntasks=64
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu=200MB
#SBATCH --output=EstimNetDirected-ifd_inouye_pyke_pollinators_altk4cycles-%j.out
#SBATCH --error=EstimNetDirected-ifd_inouye_pyke_pollinators_altk4cycles-%j.err

echo -n "started at: "; date

#ROOT=${HOME}/EstimNetDirected
ROOT=../../../

module purge
module load gcc/11.3.0    # need this for 'toolchain' hierarchical modules
module load openmpi/4.1.4 # module version numbers are required on OzStar


cp -p ${ROOT}/src/EstimNetDirected_mpi_arrays .

time mpirun ./EstimNetDirected_mpi_arrays config_inouye_pyke_pollinators_ifd_altk4cycles.txt

module load r/4.2.1 # module version numbers are required on OzStar


time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_ifd_inouye_pyke_pollinators_altk4cycles dzA_ifd_inouye_pyke_pollinators_altk4cycles |tee estimation_ifd_inouye_pyke_pollinators_altk4cycles.txt
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_ifd_inouye_pyke_pollinators_altk4cycles dzA_ifd_inouye_pyke_pollinators_altk4cycles

echo -n "ended at: "; date

