#!/bin/bash

#SBATCH --job-name="EstimNetDirected_mpi-ifd_inouye_pyke_pollinators"
#SBATCH --ntasks=20
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu=600MB
#SBATCH --output=EstimNetDirected-ifd_inouye_pyke_pollinators-%j.out
#SBATCH --error=EstimNetDirected-ifd_inouye_pyke_pollinators-%j.err

echo -n "started at: "; date

#ROOT=${HOME}/EstimNetDirected
ROOT=../../../

module purge
module load gcc/11.3.0    # need this for 'toolchain' hierarchical modules
module load openmpi/4.1.4 # module version numbers are required on OzStar


cp -p ${ROOT}/src/EstimNetDirected_mpi_arrays .

time mpirun ./EstimNetDirected_mpi_arrays config_inouye_pyke_pollinators_ifd.txt

module load r/4.2.1 # module version numbers are required on OzStar


time Rscript ${ROOT}/scripts/computeEstimNetDirectedCovariance.R theta_ifd_inouye_pyke_pollinators dzA_ifd_inouye_pyke_pollinators |tee estimation_ifd_inouye_pyke_pollinators.txt
time Rscript ${ROOT}/scripts/plotEstimNetDirectedResults.R theta_ifd_inouye_pyke_pollinators dzA_ifd_inouye_pyke_pollinators

echo -n "ended at: "; date

