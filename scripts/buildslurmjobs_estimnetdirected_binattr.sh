#!/bin/sh
#
# File:    buildslurmjobs_estimnetdirected_binattr.sh
# Author:  Alex Stivala
# Created: December 2014
#
#
# buildslurmjobs_estimnetdirected.sh - build directory tree for set of estimations 
#                                      with binary attributes
#
# Usage: buildslurmjobs_estimnetdirected.sh  out_root_dir binattr_file in_networks_files_basename
#
# Build directory heirarchy at out_root_dir with settings file and slurm
# script in subdirs for each network in sample statistics (PNet simulation
# outpt files) pattern in_networks_files
#
# Uses templatelate files and script in its own sripts directory
# The following variables are replaced:
#   @JOBNAME
#   @SETTINGS
#   @ROOT
#   @NETWORK
#   @BINARYATTRIBUTES
# in the following template files:
#   config.template
#   estimnetdirected_slurm_script.template
#

SCRIPTDIR=`dirname $0`
PBS_TEMPLATE=${SCRIPTDIR}/estimnetdirected_slurm_script.template
SETTINGS_TEMPLATE=${SCRIPTDIR}/config_binattr.template

root=`dirname ${SCRIPTDIR}`

if [ $# -ne 3 ]; then
  echo "usage: $0 out_root_dir binattr_file in_networks_files_Basename" >&2
  exit 1
fi
out_root_dir=$1
binattr_file=$2
infile_basename=$3

if [ ! -d $out_root_dir ]; then
  mkdir $out_root_dir
fi

for samplefile in ${infile_basename}*.txt
do
  samplenum=`echo ${samplefile} | sed "s!${infile_basename}!!g" | sed s'/[.]txt//g'`
  outdir=${out_root_dir}/sample${samplenum}
  mkdir ${outdir}
  networkfile=${outdir}/arclist.txt
  fgrep -i '*arcs' ${samplefile} >/dev/null 2>&1
  if [ $? -eq 0 ]; then
    # already Pajek arc list, just copy
    cp ${samplefile} ${networkfile}
  else
    # convert from matrix to arc list
    $SCRIPTDIR/extractnetwork.sh ${samplefile} | Rscript ${SCRIPTDIR}/convertMatrixToArclist.R > ${networkfile}
  fi
  binattrfilename=binaryAttributes.txt
  cat ${binattr_file} > ${outdir}/${binattrfilename}
  settingfile=${outdir}/config.txt
  networkfilename=`basename ${networkfile}`
  settingfilename=`basename ${settingfile}`
  cat $SETTINGS_TEMPLATE | sed "s!@NETWORK!${networkfilename}!g" | sed "s!@BINARYATTRIBUTES!${binattrfilename}!g"  > ${settingfile}
  jobfilename=${outdir}/estimnetdirected_mpi_slurm_script.sh
  cat $PBS_TEMPLATE | sed "s!@ROOT!${root}!g" |  sed "s!@JOBNAME!estimnetdirected_${samplenum}!g" | sed "s!@SETTINGS!${settingfilename}!g" > ${jobfilename}
done

