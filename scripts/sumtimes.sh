#!/bin/sh
#
# File:    sumtimes.sh
# Author:  Alex Stivala
# Created: August 2008
#
#
#
# sumtimes.sh - add up all the elapsed real time from time output from mpi runs
#
# Usage: sumtimes.sh [-m|-h|-s] file_list
#
# -m : output only minutes and seconds, no hours
# -h : only output hours and minutes, no seconds
# -s : only output seconds
#
# file_list list of files containing in them time output
# assume time is in bash builtin time format e.g.
#
#real    48m8.320s
#user    34m50.430s
#sys     7m20.210s
#
# Allows multiple times in each file, adding them all up.
# Output, sum of all the elpased times, is to stdout.
#
#

if [ $# -lt 1 ]; then
    echo "Usage: $0 [-m|-h] file list" >&2
    exit 1
fi

use_hours=1
use_seconds=1
only_seconds=0

if [ $# -ge 1 ] ; then
    if [ `expr substr $1 1 1` = "-" ]; then
        if [ "$1" = "-m" ]; then
            use_hours=0
        elif [ "$1" = "-h" ]; then
            use_seconds=0
        elif [ "$1" = "-s" ]; then
            only_seconds=1
        else
            echo "Usage: $0 [-m|-h|-s] file list" >&2
            exit 1
        fi
        shift 1
    fi
fi

total_seconds=0
builtinformat=0
for errfile in $*
do
   for elapsed in `grep --text '^real' ${errfile} | awk '{print $2}'`
   do
      mindex=`expr index ${elapsed} 'm'`
      mindex=`expr $mindex - 1`
      mins=`expr substr ${elapsed} 1 ${mindex}`
      secindex=`expr $mindex + 2`
      sindex=`expr index ${elapsed} 's'`
      seclen=`expr ${sindex} - ${secindex}`
      secs=`expr substr ${elapsed} ${secindex} ${seclen}`
      total_seconds=`echo "$total_seconds + $mins * 60 + $secs" | bc -l`
   done
done

total_seconds=`printf "%.0f" $total_seconds`

if [ $only_seconds -eq 1 ]; then
    printf '%d' ${total_seconds}
elif [ $use_hours -eq 1 ]; then
    hours=`expr $total_seconds / 3600`
    mins=`expr $total_seconds - $hours \* 3600`
    mins=`expr $mins / 60`
    rsecs=`expr $total_seconds - $hours \* 3600`
    rsecs=`expr $rsecs - $mins \* 60`
    if [ $use_seconds -eq 1 ]; then
        printf '%d h %02d m %02d s' ${hours} ${mins} ${rsecs}
    else
        if [ $rsecs -ge 30 ]; then
            mins=`expr $mins + 1`
        fi
        printf '%d h %02d m' ${hours} ${mins}
    fi
else
    mins=`expr $total_seconds / 60`
    rsecs=`expr $total_seconds - $mins \* 60`
    printf '%02d m %02d s' ${mins} ${rsecs}
fi
