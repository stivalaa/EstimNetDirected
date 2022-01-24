#!/bin/sh

# Regression test for two-path matrix construction/update and change
# statistics.

rc=0

BASELINE=polblogs_test_results_baseline.txt

OUTPUT=polblogs_test_results.out
DIFFILE=polblogs_test_results.diff

echo "1. no two-path lookup"

time ./testChangeStatsDirected ../../pythonDemo/polblogs/polblogs_arclist.txt  polblogs_nodepairs.txt | fgrep -v nnz | fgrep -v DEBUG  > ${OUTPUT}

BASELINE_NO2PATHTABLES=polblogs_test_results_baseline_no2pathtables.txt
cat ${BASELINE} | sed '/^[a-z]*2p sum = .*/d;s/, [a-z]*2p sum = .*$//'  > ${BASELINE_NO2PATHTABLES}
diff ${BASELINE_NO2PATHTABLES} ${OUTPUT} > ${DIFFILE}

if [ $? -eq 0 ]; then
  echo
  echo "PASSED"
else 
  echo
  echo "**** FAILED ****"
  echo "diff results are in ${DIFFILE}"
  rc=1
fi

OUTPUT=polblogs_test_results_array.out
DIFFILE=polblogs_test_results_array.diff

echo "2. two-path matrices"

time ./testChangeStatsDirected_array ../../pythonDemo/polblogs/polblogs_arclist.txt  polblogs_nodepairs.txt | fgrep -v nnz | fgrep -v DEBUG  > ${OUTPUT}

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

echo "3. two-path hash tables"

OUTPUT=polblogs_test_results_hash.out
DIFFILE=polblogs_test_results_hash.diff

time ./testChangeStatsDirected_hash ../../pythonDemo/polblogs/polblogs_arclist.txt  polblogs_nodepairs.txt | fgrep -v nnz | fgrep -v DEBUG  > ${OUTPUT}

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
