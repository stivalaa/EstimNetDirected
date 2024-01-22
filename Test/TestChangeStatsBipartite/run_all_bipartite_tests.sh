#!/bin/sh
#
# File:    run_all_bipartite_tests.sh
# Author:  Alex Stivala
# Created: May 2022 
#
# Run all the regression  / unit tests for bipartite graphs.
#

./run_test_bpnet_A6000_B750.sh
./run_test_robertson_pollinators.sh
./run_test_inouye_pyke_pollinators.sh
./run_test_stats_sum_change_stats_bipartite.sh
