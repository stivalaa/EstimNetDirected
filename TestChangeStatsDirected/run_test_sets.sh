#!/bin/sh
#
# File:    run_test_sets.sh
# Author:  Alex Stivala
# Created: July 2019
#
#
# run_test_sets.sh - regression test for set parsing and similarity 
#
# Usage: run_test_sets.sh 
#

rc=0

INPUT=../examples/setAttributes_n1000.txt
TMP_INPUT=set_test_input.txt
OUTPUT=set_test_results.out
CHECK_OUTPUT=set_test_results_check.out


cat ${INPUT} | cut -d' ' -f1 > ${TMP_INPUT}

echo "Running tests on set parsing and set similarity functions..."

time ./testSetFunctions ${TMP_INPUT} > ${OUTPUT}
python3 ./checkSetFunctions.py ${TMP_INPUT} ${OUTPUT} > ${CHECK_OUTPUT}

if [ $? -eq 0 ]; then
    echo
    echo "PASSSED"
else 
  echo
  echo "**** FAILED ****"
  echo "check results  in ${CHECK_OUTPUT}"
  rc=1
fi

exit $rc
