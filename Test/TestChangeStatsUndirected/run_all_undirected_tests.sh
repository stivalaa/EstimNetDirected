#!/bin/sh
#
# File:    run_all_undirected_test.sh
# Author:  Alex Stivala
# Created: January 2024
#
# Run all the regression  / unit tests for undirected graphs.
#
fail=0
./run_test_netscience.sh || fail=1
./run_test_stats_sum_change_stats_undirected.sh || fail=1
./run_test_diff_stats_change_stats_undirected.sh || fail=1
if [ $fail -ne 0 ]; then
    echo "***** A test in $0 FAILED *****"
fi
exit $fail

