#!/bin/sh
#
# File:    reorderEstimnetTextable.sh
# Author:  Alex Stivala
# Created: December 2017
#
#
# Read output of estimnetdirectedEstimation2textableMultiModels.sh and permute the table
# rows according to specified order, in order to present the rows in an order we determine
# as best for presentation.
# 
# Usage: reorderEstimnetTextable.sh <intable> <permutation>
#
# where <intable> is the input tex file and
#       <permutation> is a permutation of 0..n-1 e.g. 2 0 1
#
# E.g.:
#   reorderEstimnetTextable.sh  ecoli_estimations.tex  1 5 4 2 3 0 9 10 7 6
#
# Output is to stdout
#
# Uses various GNU utils options on echo, etc.


# https://stackoverflow.com/questions/15639888/reorder-lines-of-file-by-given-sequence
# Usage: schwartzianTransform "A.txt" 2 0 1
function schwartzianTransform {
    local file="$1"
    shift
    local sequence="$@"
    echo -n "$sequence" | sed 's/[^[:digit:]][^[:digit:]]*/\
/g' | paste -d ' ' - "$file" | sort -n | sed 's/^[[:digit:]]* //'
}

if [ $# -lt 2 ]; then
  echo "Usage: $0 <infile> <permutation>" >&2
  exit 1
fi
infile=$1
shift 1
permutation=$*


tmpfile=`mktemp`

nfields=`echo ${permutation} | wc -w`
nuniq=`echo ${permutation} | tr ' ' '\n' | sort -n | uniq | wc -w`

if [ $nuniq != $nfields ]; then
  echo Not a valid permutation: $nfields fields but $nuniq unique fields >&2
  exit 1
fi

# assume format of input latex table from estimnetdirectedEstimation2textableMultiModels.sh
# 7 header rows and 2 trailer rows
headrows=7
tailrows=2
head -n${headrows} ${infile}
cat ${infile} | tail -n+`expr ${headrows} + 1` |head -n-${tailrows} > ${tmpfile}
nrows=`cat ${tmpfile} | wc -l`
if [ $nrows -ne $nfields ]; then
  echo Permutation has $nfields fields but table has $nrows rows >&2
  rm ${tmpfile}
  exit 1
fi
schwartzianTransform ${tmpfile} ${permutation} 
tail -n${tailrows} ${infile}

rm ${tmpfile}

