#!/usr/bin/env python2
##############################################################################
#
# genRandomNodeNumsBipartite.py - generate random node pairs for two-mode graph
#
# Generate random distinct node number pairs i,j 
# with i in range 0 .. Np-1 and j in range 0 .. N-Np-1
# where Np is number of 'type P' (rows) nodes and N is total number of nodes
#
# File:    genRandomNodeNumsBipartite.py
# Author:  Alex Stivala
# Created: April 2022
#
##############################################################################

"""
genRandomNodeNumsBipartite.py - generate random node pairs for two-mode graph

Generate random distinct node number pairs i,j 
with i in range 0 .. Np-1 and j in range 0 .. N-Np-1
where Np is number of 'type P' (rows) nodes and N is total number of nodes

For use with two-mode networks in Pajek format:

   the first lines should be e.g.
   *vertices 36 10
   first number is total number of nodes
   second number is number of type P nodes ['people'].
   the rest are type A ['affiliation'] - conventionally in the affiliation
   matrix the rows are people (P) and the columns affiliations (A).
   They must be numbered 1 ... N where N = num_vP + num_vA
   so nodes 1 .. num_vP are type P and num_vP+1 .. N are type A

Usage:
    genRandomNodeNums.py num_pairs N Np
    
    num_pairs is the number of pairs to generate
    N is total number of nodes
    Np is number of 'type P' (rows) nodes
   
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
    sys.stderr.write("usage: " + progname + " num_pairs N Np\n")
    sys.stderr.write("   num_pairs is the number of node pairs to generate\n")
    sys.stderr.write("   N is the total number of nodes\n")
    sys.stderr.write("   Np is the number of type P nodes \n")
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

    if len(args) != 3:
        usage(sys.argv[0])

    num_pairs = int(args[0])
    N = int(args[1])
    Np = int(args[2])

    if Np > N:
      sys.stderr.write("ERROR: Must have Np <= N\n")
      sys.exit(1)
    
    for i in xrange(num_pairs):
      nodepair_i = random.randint(0, Np-1)
      nodepair_j = random.randint(0, N-Np-1)
      sys.stdout.write("%d %d\n" % (nodepair_i, nodepair_j))

    
if __name__ == "__main__":
    main()


