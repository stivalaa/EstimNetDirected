#!/usr/bin/env python3
##############################################################################
#
# genRandomNodeNums.py - generate random continuous attributes
#
# Generate random distinct node number pairs each in range 0..N-1
#
# File:    genRandomNodeNums.py
# Author:  Alex Stivala
# Created: April 2018
#
##############################################################################

"""
Generate distinct pairs of integers in range 0..N-1

Usage:
    genRandomNodeNums.py num_pairs N
    
    num_pairs is the number of pairs to generate
    N is number of nodes (so generate numbers in [0,N-1])
   
Output is to stdout, one space-separated pair per line.
"""

import sys
import getopt
import random


#-----------------------------------------------------------------------------
#
# Functions
#
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#
# main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname + " num_pairs N\n")
    sys.stderr.write("   num_pairs is the number of node pairs to generate\n")
    sys.stderr.write("   N is the number of nodes (so generates in [0..N-1])\n")
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

    num_pairs = int(args[0])
    N = int(args[1])
    
    for i in range(num_pairs):
      nodepair = random.sample(range(N), 2)
      sys.stdout.write("%d %d\n" % (nodepair[0], nodepair[1]))

    
if __name__ == "__main__":
    main()


