#!/bin/bash

#SBATCH --job-name="Py_EstimNetSimple"
#SBATCH --ntasks=1
#SBATCH --time=0-8:00:00

echo -n "started at: "; date

module load numpy/1.8.2-vlsci_intel-2015.08.25-Python-2.7.9
module load R/3.2.1-vlsci_intel-2015.08.25

time python ./runExample.py

time Rscript plotEstimNetSimpleDemoResults.R theta_values_sample_statistics_n500_directed_binattr_sim420000000.txt dzA_values_sample_statistics_n500_directed_binattr_sim420000000.txt

echo -n "ended at: "; date

