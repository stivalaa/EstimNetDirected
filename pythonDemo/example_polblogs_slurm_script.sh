#!/bin/bash

#SBATCH --job-name="Py_EstimNetSimple"
#SBATCH --ntasks=1
#SBATCH --time=0-48:00:00

echo -n "started at: "; date

module load python/3.8.5

time python3 ./runExample_polblogs.py

module unload python # on cluster module load r will not work if this is not done
module load r

time Rscript plotEstimNetSimpleDemoResults.R theta_values_polblogs_arclist.txt dzA_values_polblogs_arclist.txt

echo -n "ended at: "; date

