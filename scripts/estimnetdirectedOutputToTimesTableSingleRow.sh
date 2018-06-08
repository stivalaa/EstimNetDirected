#!/bin/sh
#
# File:    estimnetdirectedOutputToTimesTableSingleRow.sh
# Author:  Alex Stivala
# Created: December 2017
#
#
# Read N, density, times from stdout (slurm job file output) of 
# EstimNetDirected and format into LaTeX table.
# 
#
# Usage: estimnetdirectedOutputToTimesTableSingleRow.sh networkName ModelName slurm-file.out
#
# E.g.:
#    estimnetdirectedOutputToTimesTableSingleRow.sh "Hospital" "Model 1" ../hospital/slurm-3300854.out
#
# Output is to stdout.
#
# Uses various GNU utils options on  echo &  etc.


if [ $# -ne 3 ] ; then
  echo "Usage: $0 NetworkName ModelName slurm-file" >&2
  exit 1
fi
networkName=$1
modelName=$2
slurmfile=$3

nodecount=`fgrep -w density $slurmfile | awk '{print $3}'`
density=`fgrep -w density $slurmfile | awk '{print $9}' | tr -d ')'`
elapsed=`sumtimes.sh -m $slurmfile`

printf "%s & %d & %.5f & %s & %s \\\\\\ \n" "${networkName}" $nodecount $density "${modelName}" "${elapsed}"

