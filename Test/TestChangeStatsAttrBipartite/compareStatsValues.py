#!/usr/bin/env python3
##############################################################################
#
# compareStatsValues.py - compare values of stats computed by two methods
#
# Read colums of stats produced by statnet and columns of values produced
# by EstimNetDirected stats values test program and compare them for
# floating point equality.
#
# File:    compareStatsValues.py
# Author:  Alex Stivala
# Created: December 2023
#
##############################################################################

"""
  compareStatsValues.py - compare values of stats computed by two methods

  Read colums of stats produced by statnet and columns of values produced
  by EstimNetDirected stats values test program and compare them for
  floating point equality.

  Usage:
    compareStatsValues.py baselineoutput.txt testoutput.txt


"""

import sys,getopt
from math import isclose
    
#-----------------------------------------------------------------------------
#
# main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname + " baselineoutput.txt testoutput.txt\n")
    sys.exit(1)

def main():
    """
    See usage message in module header block
    """
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        usage(sys.argv[0])

    if len(args) != 2:
        usage(sys.argv[0])

    exitcode = 0

    baselinetext = args[0]
    testtext = args[1]

    

    # Note in the following,
    #  map(list, zip(*[row.split() for row in open(filename).readlines()]))
    # reads the data and transposes it so we have a list of columns
    # not a list of rows, which then makes it easy to convert to
    # the dict where key is column header and value is list of values
    # (converted to the appropriate data type from sting with float_or_na()
    # etc.)
    # https://stackoverflow.com/questions/6473679/transpose-list-of-lists#
    

    ## For baseline (statnet) ouitput, filter out the names from summary
    ## output (that all start with 'b', e.g. "b1nodematch.catattrA" etc.
    ## so we just get the statistic values (Which therefore must be
    ## in same order as those output by test program)
    baseline_values = list(list(map(float, col)) for col in map(list, list(zip(*[row.split() for row in [s for s in open(baselinetext).readlines() if not s.startswith("b")]]))))

    test_values = list(list(map(float, col)) for col in map(list, list(zip(*[row.split() for row in open(testtext).readlines()]))))
    
    #print( baseline_values)
    #print(test_values)

    assert len(baseline_values) == len(test_values)

    for i in range(len(baseline_values)):
        isequal = map(lambda x,y : isclose(x, y, rel_tol=0.0001), baseline_values[i], test_values[i])
        if not all(isequal):
            print("FAILED at column", i, " row(s) ", [j for j,x in enumerate(isequal) if x], ":")
            print(baseline_values[i])
            print(test_values[i])
            sys.exit(1)
    print("PASSED")
    
    
if __name__ == "__main__":
    main()


