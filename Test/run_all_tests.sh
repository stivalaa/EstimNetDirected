#!/bin/sh
#
# File:    run_all_tests.sh
# Author:  Alex Stivala
# Created: July 2019
#
# Run all the regression  / unit tests.
#
fail=0
(cd TestChangeStatsDirected && ./run_test_polblogs.sh) || fail=1
(cd TestChangeStatsUndirected && ./run_all_undirected_tests.sh) || fail=1
(cd TestSetFunctions && ./run_test_sets.sh) || fail=1
(cd TestChangeStatsBipartite && ./run_all_bipartite_tests.sh) || fail=1
(cd TestChangeStatsAttrBipartite && ./run_all_attr_bipartite_tests.sh) || fail=1
##(cd TestRandom && ./run_test_random.sh) || fail=1
if [ $fail -ne 0 ]; then
    echo "***** A test FAILED *****"
fi
exit $fail
