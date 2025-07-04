#!/bin/sh
#
# File:    estimnetdirectedEstimation2tableMultiModels.sh
# Author:  Alex Stivala
# Created: August 2021
#
# (Based on estimnetdirectedEStimation2textableMultiModels.sh)
#
# Read output of computeEstimNetDirectedCovariance.R with the estimate,
# estimated std. error and t-ratio computed from EstimNetDirected results
# and write plain text table (whitepace spearte file) for convenient reading in
# R, for multiple different models
# 
# Usage: estimnetdirectedEstimation2tableMultiModels.sh estimationoutputfile_model1 estimationoutputfile_model2 ...
#
# E.g.:
#   estimnetdirectedEstimation2tableMultiModels.sh  estimation.out model2/estimation.out
#
# Output is to stdout
#
# Uses various GNU utils options on echo, etc.


if [ $# -lt 1 ]; then
    echo "usage: $0 estimation1.out estimation2.out ..." >&2
    exit 1
fi

num_models=`expr $#`

estimnet_tmpfile=`mktemp`

echo "# Generated by: $0 $*"
echo "# At: " `date`
echo "# On: " `uname -a`

echo model parameter estimate stderr tratio

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
model=1
for estimationresults in $*
do
    # https://unix.stackexchange.com/questions/78472/print-lines-between-start-and-end-using-sed
    cat ${estimationresults} | sed -n -e '/^Pooled/,${//!p}'  | tr -d '*' | fgrep -vw AcceptanceRate | awk '{print $1,$2,$4,$5}'  |  tr ' ' '\t' | sed "s/^/${model}\t/" >> ${estimnet_tmpfile}
    model=`expr $model + 1`
done


effectlist=`cat ${estimnet_tmpfile} | fgrep -vw TotalRuns | fgrep -vw ConvergedRuns | awk '{print $2}' | sort | uniq`

for model in `seq 1 $num_models`
do
    for effect in ${effectlist} ConvergedRuns TotalRuns
    do
        echo -n "${model} ${effect} " 
        if [ ${effect} = "ConvergedRuns" -o ${effect} = "TotalRuns" ]; then
               runs=`grep -w ${effect} ${estimnet_tmpfile} | awk -vmodel=$model '$1 = model {print $3}'`
               echo -n " ${runs} NA NA"
        else
            estimnet_point=`grep -w ${effect} ${estimnet_tmpfile} | awk -vmodel=$model '$1 = model {print $3}'`
            estimnet_stderr=`grep -w ${effect} ${estimnet_tmpfile} | awk -vmodel=$model '$1 = model {print $4}'`
            estimnet_tratio=`grep -w ${effect} ${estimnet_tmpfile} | awk -vmodel=$model '$1 = model {print $5}'`
            if [ "${estimnet_point}" = "" ];  then
                echo -n " NA NA NA"
            else 
                echo -n " ${estimnet_point} ${estimnet_stderr} ${estimnet_tratio}"
            fi
        fi
           echo
     done
done

rm ${estimnet_tmpfile}
