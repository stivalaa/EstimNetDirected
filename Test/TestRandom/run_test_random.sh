#!/bin/bash
#
# File:    run_test_sets.sh
# Author:  Alex Stivala
# Created: June 2022
#
#
# run_test_random.sh - regression test for random integers
#
# Usage: run_test_sets.sh 
#

set -o pipefail

OUTPUT=random_test_results.out
SAMPLE_SIZE=100000000

echo "Running tests on random integer generation..."
echo -n > ${OUTPUT}
for n in 1000 10000 1000000 10000000
do
  time ./testRandom ${n} ${SAMPLE_SIZE} | Rscript checkUniform.R ${n} | tee -a ${OUTPUT}
  if [ $? -ne 0 ]; then
    echo "**** FAILED ****"
    exit 1
  fi
done
echo "PASSED"
exit 0

