#!/usr/bin/env python
##############################################################################
#
# checkSetFunctions.py - verify set similarity (Jaccard) results
#
# Parses sets in EstimNEtDirected format and output file from
# run_test_sets.sh and checks that results are correct.
#
# File:    checkSetFunctions.py
# Author:  Alex Stivala
# Created: July 2019
#
##############################################################################

"""
Verify set similarity results.

Usage:
    checkSetFunctions.py setattributes testoutput
    
    setattributes is input file (Created by run_test_sets.sh)
    testoutput is output file from run_test_sets.sh runnig ntestSetFunctions
    
Output is to stdout 
"""

import sys
import getopt
import random


#-----------------------------------------------------------------------------
#
# Functions
#
#-----------------------------------------------------------------------------

def jaccard_index(a, b):
    """
    Jaccard similarity of sets a and b
    """
    intsize = len(set.intersection(a, b))
    unionsize = len(set.union(a, b))
    if unionsize == 0:
        return 1
    else:
        return float(intsize) / float(unionsize)
    
#-----------------------------------------------------------------------------
#
# main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname + " setinput testoutput\n")
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

    setinput = args[0]
    testoutput = args[1]

    setlist = [set() if s.rstrip() == 'none' else set(s.rstrip().split(',')) for s in open(setinput).readlines()[1:]]

    for line in open(testoutput).readlines():
        (i, j, sim) = line.split()
        i = int(i)
        j = int(j)
        sim = float(sim)
        vsim = jaccard_index(setlist[i], setlist[j])
        if (abs(sim - vsim) > 1e-05):
            print 'FAIL',i,j,sim,vsim
            exitcode = -1

    sys.exit(exitcode)
    
    
if __name__ == "__main__":
    main()


