#!/bin/sh

# Regression test for two-path matrix construction/update and change
# statistics.

rc=0

BASELINE=polblogs_test_results_baseline.txt
OUTPUT=polblogs_test_results.out
DIFFILE=polblogs_test_results.diff

echo "1. two-path matrices"

time ./testChangeStatsDirected ../pythonDemo/polblogs/polblogs_arclist.txt  polblogs_nodepairs.txt | fgrep -v nnz | fgrep -v DEBUG  > ${OUTPUT}

diff ${BASELINE} ${OUTPUT} > ${DIFFILE}

if [ $? -eq 0 ]; then
  echo
  echo "PASSED"
else 
  echo
  echo "**** FAILED ****"
  echo "diff results are in ${DIFFILE}"
  rc=1
fi

echo "2. two-path hash tables"

OUTPUT=polblogs_test_results_hash.out
DIFFILE=polblogs_test_results_hash.diff

time ./testChangeStatsDirected_hash ../pythonDemo/polblogs/polblogs_arclist.txt  polblogs_nodepairs.txt | fgrep -v nnz | fgrep -v DEBUG  > ${OUTPUT}

diff ${BASELINE} ${OUTPUT} > ${DIFFILE}

if [ $? -eq 0 ]; then
  echo
  echo "PASSED"
else 
  echo
  echo "**** FAILED ****"
  echo "diff results are in ${DIFFILE}"
  rc=2
fi

exit $rc
