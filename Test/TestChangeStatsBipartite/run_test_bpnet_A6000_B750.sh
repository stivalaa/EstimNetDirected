#!/bin/sh

# Regression test for two-path matrix construction/update and change
# statistics.

rc=0

BASELINE=bpnet_A6000_B750_test_results_baseline.txt

OUTPUT=bpnet_A6000_B750_test_results.out
DIFFILE=bpnet_A6000_B750_test_results.diff

echo "1. no two-path lookup"

time ./testChangeStatsBipartite ../../examples/bipartite/simulated/bpnet_A6000_B750_sparse_sim100000000.net  bpnet_A6000_B750_nodepairs.txt | fgrep -v nnz | fgrep -v DEBUG  > ${OUTPUT}

BASELINE_NO2PATHTABLES=bpnet_A6000_B750_test_results_baseline_no2pathtables.txt
cat ${BASELINE} | sed '/^[A-Za-z]*2p sum = .*/d;s/, [A-Za-z]*2p sum = .*$//'  > ${BASELINE_NO2PATHTABLES}
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

OUTPUT=bpnet_A6000_B750_test_results_array.out
DIFFILE=bpnet_A6000_B750_test_results_array.diff

echo "2. two-path matrices"

time ./testChangeStatsBipartite_array ../../examples/bipartite/simulated/bpnet_A6000_B750_sparse_sim100000000.net  bpnet_A6000_B750_nodepairs.txt | fgrep -v nnz | fgrep -v DEBUG  > ${OUTPUT}

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

OUTPUT=bpnet_A6000_B750_test_results_hash.out
DIFFILE=bpnet_A6000_B750_test_results_hash.diff

time ./testChangeStatsBipartite_hash ../../examples/bipartite/simulated/bpnet_A6000_B750_sparse_sim100000000.net  bpnet_A6000_B750_nodepairs.txt | fgrep -v nnz | fgrep -v DEBUG  > ${OUTPUT}

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
