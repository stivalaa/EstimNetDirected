#!/bin/bash

#SBATCH --job-name="Py_EstimNetSimple"
#SBATCH --ntasks=1
#SBATCH --time=0-24:00:00

echo -n "started at: "; date

module load numpy/1.8.2-vlsci_intel-2015.08.25-Python-2.7.9
module load R/3.2.1-vlsci_intel-2015.08.25

time python ./runExample_polblogs.py

time Rscript plotEstimNetSimpleDemoResults.R theta_values_polblogs_arclist.txt dzA_values_polblogs_arclist.txt

echo -n "ended at: "; date

