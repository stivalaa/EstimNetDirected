#!/bin/bash

#SBATCH --job-name="EstimNetDirected_regression_tests"
#SBATCH --ntasks=1
#SBATCH --time=0-02:00:00
#SBATCH --mem-per-cpu=2GB
#SBATCH --output=EstimNetDirected_regression_tests-%j.out
#SBATCH --error=EstimNetDirected_regression_tests-%j.err

echo -n "started at: "; date

module purge
module load gcc/11.3.0     # need this for 'toolchain' hierarchical modules
module load openmpi/4.1.4 # module version numbers are required on OzStar
module load r/4.2.1

time ./run_all_tests.sh

echo -n "ended at: "; date

