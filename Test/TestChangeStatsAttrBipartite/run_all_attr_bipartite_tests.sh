#!/bin/sh
#
# File:    run_all_attr_bipartite_tests.sh
# Author:  Alex Stivala
# Created: January 2024
#
# Run all the regression  / unit tests for bipartite graphs with attributes
#
fail=0
./run_test_bpnet_A12000_B4000_attr.sh || fail=1
./run_test_b1nodematch_bpnet_A12000_B4000_attr.sh || fail=1
./run_test_stats_sum_change_stats_attr_bipartite.sh || fail=1
if [ $fail -ne 0 ]; then
    echo "***** A test in $0 FAILED *****"
fi
exit $fail

