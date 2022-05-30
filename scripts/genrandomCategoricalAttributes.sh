#!/bin/sh
#
# genrandomCategoricalAttributes.sh
#
# Generate random categorical attributes, uniformly distributed.
#
# Usage:
#    genrandomCategoricalAttributes N q
#
#      N is the number of attributes to generate
#      q is the number of different categories (generated as 0..q-1)
#
# Output is to stdout in PNet format for attributes; first line is header
# 
# ADS 4June2014
#
# Uses GNU utils seq, shuf etc.
#
if [ $# -ne 2 ]; then
  echo "Usage: $0 N q" >&2
  exit 1
fi
N=$1
q=$2

echo "categoricalAttribute"
for i in `seq $N`
do
  echo `expr $i \% $q`
done | shuf | sed 's/$/\r/'
# NB must have unix2dos (now using sed to do it as unix2dos not present on 
# some systems) otherwise PNet crashes jvm with no error message


