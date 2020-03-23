#!/bin/sh
#
# File:    estimnetdirectedEstimation2simulationConfig.sh
# Author:  Alex Stivala
# Created: November 2019
#
#
# Read output of computeEstimNetDirectedCovariance.R with the estimate
# computed from EstimNetDirected results
# and build config script for SimulateERGM to simulate networks from
# the estimated parmeter values. Also needs to parse some filenames
# from the configuration file used for the estimation, to get arc list
# file (for number of nodes), output file names (obseved network statistics)
# etc.
# 
# Usage: estimnetdirectedEstimation2simulationConfig.sh estimation_config_file estimationoutputfile statsoutputfile simNetFilePrefix
#
#  estmiation_config_file is the config file that generated the
#      estimationoutputfile
#  statsoutputfile is the file to write the simulated network stats to
#    (this is written in the output config file)
#  simNetFilePrefix is prefix for simulated Pajek .net files
#    (this is written in the output config file)
#
# E.g.:
#   estimnetdirectedEstimation2simulationConfig.sh config_example.txt estimation.out stats_estimation.out gof_estimation
#
# Output is to stdout
#
# Uses various GNU utils options on echo, etc.



if [ $# -ne 4 ]; then
    echo "usage: $0 estimation_config.txt estimation.out statsoutputfilename simNetFilePrefix" >&2
    exit 1
fi

estimationconfig=$1
estimationresults=$2
statsFile=$3
simNetFilePrefix=$4

estimnet_tmpfile=`mktemp`
estimnet_tmpfile2=`mktemp`


echo "# Generated by: $0 $*"
echo "# At: " `date`
echo "# On: " `uname -a`

arclistFile=`grep -i arclistFile ${estimationconfig} | awk -F= '{print $2}'`
observedStatsFilePrefix=`grep -i observedStatsFilePrefix ${estimationconfig} | awk -F= '{print $2}'`

echo "# arclistFile = ${arclistFile}"
echo "# observedStatsFilePrefix = ${observedStatsFilePrefix}"

numNodes=`cat ${arclistFile} | grep -i '^*Vertices'| awk '{print $2}'`

echo "numNodes = ${numNodes}"

cat <<EOF
useTNTsampler = True # use the tie-no-tie sampler
sampleSize = 100 #number of network samples to take from simulation
interval = 10000000 # interval (iterations) between samples
burnin = 100000000 # number of iterations to throw away before first sample
outputSimulatedNetworks = True
EOF

# we need the attribute files, directly from the esetimation config file
grep -i binattrFile ${estimationconfig}
grep -i catattrFile ${estimationconfig}
grep -i contattrFile ${estimationconfig}
grep -i setattrFile ${estimationconfig}

# options that need to be same in simulation as estimation
grep -i forbidreciprocity ${estimationconfig}

echo "# Filename of file to write statistics to"
echo "statsFile = ${statsFile}"
echo "# Prefix of simulated networks in Pajek .net file format"
echo "simNetFilePrefix = ${simNetFilePrefix}"

# new version has results starting at line following "Pooled" at start
# of line (pooling the individual run estimates values printed earlier) and
# 5 columns:
# Effect   estimate   sd(theta)   est.std.err  t.ratio
# (and maybe *) plus
# TotalRuns and ConvergedRuns e.g.:
#Diff_completion_percentage -0.002270358 0.005812427 0.01295886 0.021386
#TotalRuns 2
#ConvergedRuns 2
# (see computeEstimNetDirectedCovariance.R)
# https://unix.stackexchange.com/questions/78472/print-lines-between-start-and-end-using-sed
cat ${estimationresults} | sed -n -e '/^Pooled/,/^TotalRuns/{//!p}'  | tr -d '*' | fgrep -vw AcceptanceRate | fgrep -vw TotalRuns | fgrep -vw ConvergedRuns | awk '{print $1,$2,$4,$5}'  |  tr ' ' '\t' >> ${estimnet_tmpfile}

effectlist=`cat ${estimnet_tmpfile} |  awk '{print $1}' | sort | uniq`

for effect in ${effectlist}
do
  estimnet_point=`grep -w ${effect} ${estimnet_tmpfile} | awk '{print $2}'`
  echo $effect = $estimnet_point
done > ${estimnet_tmpfile2}

# Note this separation of the structural from attribute effects
# depends on assuming that only the attribute effects have an underscore
# in the names  (separating the parameter from the attribute name) so
# will break if any of the attributes have an underscore in them... 

echo 'structParams = {'
cat ${estimnet_tmpfile2} | fgrep -v _ | sed 's/$/,/' | tr -d '\n' | sed 's/,$/}/' | sed 's/,/,\n/g'
echo
echo
# only do attrParams if any in estimation
fgrep -q _ ${estimnet_tmpfile2}
if [ $? -eq 0 ]; then
  echo 'attrParams = {'
  # convert e.g. "MatchingReciprocity_value = -1.712176"
  # to "MatchingReciprocity(value = -1.712176)"
  cat ${estimnet_tmpfile2} | fgrep _ | sed 's/\([a-zA-Z]*\)_\([a-zA-Z0-9_]*\) = \([0-9.e-]*\)/\1(\2 = \3)/g' | sed 's/$/,/' | tr -d '\n' | sed 's/,$/}/' | sed 's/,/,\n/g'
  echo
fi

rm ${estimnet_tmpfile} ${estimnet_tmpfile2}