#!/usr/bin/env python
################################################################################
#
# File:    bipartite2adjmatrix.py
# Author:  Alex Stivala
# Created: July 2012
#
################################################################################

"""
Read the bipartite network rectangular matrix for an (n, m) bipartite
graph in the format of a 0-1 matrix with n rows and m columns where
Aij = 1 only for a tie between nodes i and j (the biadjacency matrix)
and convert to the standard adjacency matrix format (for general
graphs) which is a square matrix with dimension n+m.

See usage in docstring for main()

NB on Windows XP I find that something like

bipartite2adjmatrix.py < small_biparite.txt

does NOT work (get errno 9 "bad file decriptior"), you must use

python bipartite2adjmatrix.py < small_biparite.txt

instead (i.e. explicitly run the python intepreter with script as argument,
not execute the script with .py extension). See:

http://mail.python.org/pipermail/python-bugs-list/2004-August/024920.html

"""

import os
import sys
from bipartitematrix import read_bipartite_matrix,bipartite_to_adjmatrix
    
#-------------------------------------------------------------------------------
#
# Main
#
#-------------------------------------------------------------------------------

def usage(progname):
    """
    print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + " < bipartitematrix\n")
    sys.exit(1)

def main():
    """
    main for bipartite2adjmatrix.py

    Usage: bipartite2adjmatrix.py < bipartitematrix

    Read the bipartite network rectangual matrix on stdin (see module docstring)
    and write the correpsonding general (square) adjacency matrix to stdout.
    """
    
    if len(sys.argv) != 1:
        usage(os.path.basename(sys.argv[0]))

    (m, bipartite_graph) = read_bipartite_matrix(sys.stdin)
    general_graph = bipartite_to_adjmatrix(m, bipartite_graph)
    for i in xrange(len(general_graph)):
        for j in xrange(len(general_graph)):
            if general_graph[i].has_key(j):
                sys.stdout.write("1")
            else:
                sys.stdout.write("0")
            if j < len(general_graph)-1:
                sys.stdout.write(" ")
        sys.stdout.write("\n")

if __name__ == "__main__":
    main()

