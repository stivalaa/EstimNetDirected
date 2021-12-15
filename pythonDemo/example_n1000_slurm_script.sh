#!/bin/bash

#SBATCH --job-name="Py_EstimNetSimple"
#SBATCH --ntasks=1
#SBATCH --time=0-24:00:00

echo -n "started at: "; date


module load python/3.8.5

time python3 ./runExample_n1000.py

module unload python # on cluster module load r will not work if this is not done
module load r

time Rscript plotEstimNetSimpleDemoResults.R theta_values_sample_statistics_n1000_directed_binattr_sim620000000.txt dzA_values_sample_statistics_n1000_directed_binattr_sim620000000.txt

echo -n "ended at: "; date

