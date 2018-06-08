#!/bin/sh
#
# File:    estimnetdirectedOutputToTimesTableMultiModels.sh
# Author:  Alex Stivala
# Created: December 2017
#
#
# Read N, density, times from stdout (slurm job file output) of 
# EstimNetDirected and format into LaTeX table.
# 
#
# Usage: estimnetdirectedOutputToTimesTableMultiModels.sh networkName slurm-file1 slurm-file2 ...
#
# E.g.:
#    estimnetdirectedOutputToTimesTableMultiModels.sh "Hospital" ../hospital/slurm-3300854.out ../hospital/slurm-3302739.out
#
# Output is to stdout.
#
# Uses various GNU utils options on  echo &  etc.


if [ $# -lt 2 ] ; then
  echo "Usage: $0 NetworkName slurm-file1 slurm-file2 ..." >&2
  exit 1
fi
networkName=$1
shift 1
model=1
for slurmfile in $*
do
	estimnetdirectedOutputToTimesTableSingleRow.sh "${networkName}" "Model ${model}" ${slurmfile} 
  model=`expr $model + 1`
done

