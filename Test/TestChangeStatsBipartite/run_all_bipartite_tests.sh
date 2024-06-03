#!/bin/sh
#
# File:    run_all_bipartite_tests.sh
# Author:  Alex Stivala
# Created: May 2022 
#
# Run all the regression  / unit tests for bipartite graphs.
#
fail=0
./run_test_bpnet_A6000_B750.sh || fail=1
./run_test_robertson_pollinators.sh || fail=1
./run_test_inouye_pyke_pollinators.sh || fail=1
./run_test_stats_sum_change_stats_bipartite.sh || fail=1
./run_test_diff_stats_change_stats_bipartite.sh || fail=1
if [ $fail -ne 0 ]; then
    echo "***** A test in $0 FAILED *****"
fi
exit $fail

