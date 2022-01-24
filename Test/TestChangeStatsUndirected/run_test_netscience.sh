#!/bin/sh

# Regression test for two-path matrix construction/update and change
# statistics.

rc=0

BASELINE=netscience_test_results_baseline.txt

OUTPUT=netscience_test_results.out
DIFFILE=netscience_test_results.diff

echo "1. no two-path lookup"

time ./testChangeStatsUndirected netscience_edgelist.txt  netscience_nodepairs.txt | fgrep -v nnz | fgrep -v DEBUG  > ${OUTPUT}

BASELINE_NO2PATHTABLES=netscience_test_results_baseline_no2pathtables.txt
cat ${BASELINE} | sed '/^twoPath sum = .*/d;s/, twoPath sum = .*$//'  > ${BASELINE_NO2PATHTABLES}
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

OUTPUT=netscience_test_results_array.out
DIFFILE=netscience_test_results_array.diff

echo "2. two-path matrices"

time ./testChangeStatsUndirected_array netscience_edgelist.txt  netscience_nodepairs.txt | fgrep -v nnz | fgrep -v DEBUG  > ${OUTPUT}

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

OUTPUT=netscience_test_results_hash.out
DIFFILE=netscience_test_results_hash.diff

time ./testChangeStatsUndirected_hash netscience_edgelist.txt  netscience_nodepairs.txt | fgrep -v nnz | fgrep -v DEBUG  > ${OUTPUT}

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
