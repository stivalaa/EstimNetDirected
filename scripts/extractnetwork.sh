#!/bin/sh
#
# File:    extractnetwork.sh
# Author:  Alex Stivala
# Created: September 2013
#
#
# extractnetwork.sh - extract adjcacency matrix from PNet simulation output
#
# Usage: extractnetwork.sh  pnet_simulation_sample_statistics_file
#
# Read sample statistics file from PNet simulation output
# and write only the adjacency matrix from it to stdout.
#
#

if [ $# -ne 1 ]; then
  echo "usage: $0 sample_stats_file_name" >&2
  exit 1
fi
infilename=$1

num_vertices=`dos2unix < ${infilename} | grep '^\*vertices' | awk '{print $2}'`
dos2unix < ${infilename} | sed -n "/^\*matrix/,+${num_vertices}p"  | tail -n ${num_vertices}

