#!/bin/sh
#
# File:    run_all_tests.sh
# Author:  Alex Stivala
# Created: July 2019
#
# Run all the regression  / unit tests.
#

(cd TestChangeStatsDirected && ./run_test_polblogs.sh)
(cd TestChangeStatsUndirected && ./run_test_netscience.sh)
(cd TestSetFunctions && ./run_test_sets.sh)
(cd TestRandom && ./run_test_random.sh)
