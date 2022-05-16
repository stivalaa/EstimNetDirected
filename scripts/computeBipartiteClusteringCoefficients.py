#!/usr/bin/env python
##############################################################################
#
# computeBipartiteClusteringCoefficients.py - use NetworkX to compute bipartite c.c.
#
#
# File:    computeBipartiteClusteringCoefficients.py
# Author:  Alex Stivala
# Created: May 2022
#
##############################################################################

"""
Load bipartite netowrk in Pajek bipartite network format 
and compute Bipartite clustering coefficients. Using NetworkX
as igraph does not have bipartite functions (as of igraph 1.2.10) [using R]
and tnet [also R] is too slow and/or runs out of memory trying to compute them
while Python and NetworkX on the same network takes less than a second (and no noticable
memory usagE).

Using NetworkX version 2.6.3 and Python version 3.9.10

For NetworkX see https://networkx.org/

"""
import sys
import networkx as nx
from networkx.algorithms import bipartite


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname + " pajek_bipartite_net_filename\n")
    sys.exit(1)

def main():
    """
    See usage message in module header block
    """
    
    if len(sys.argv) != 2:
        usage(sys.argv[0])

    pajek_filename = sys.argv[1]

    # read_pajek() does not work for bipartite Pajek format so have to read and build
    # graph from file here
    # We assume the nodes are numbered 1..N (But renumber them from 0..N-1 for Python)

    pajek_text = open(pajek_filename).readlines()

    # first line e..g
    # *vertices 16000 12000
    # give total number of nodes and number of type A nodes (rest are type B)
    firstline = pajek_text[0].rstrip().lower().split(' ')
    if len(firstline) != 3 or firstline[0] != '*vertices':
        sys.stderr.write('expecting *vertices num_nodes num_A_nodes as first line\n')
        sys.exit(1)
    num_nodes = int(firstline[1])
    num_A_nodes = int(firstline[2])
    
    G = nx.Graph()
    G.add_nodes_from(range(0, num_A_nodes), bipartite = 0)
    G.add_nodes_from(range(num_A_nodes, num_nodes), bipartite = 1)
    # edges are after '*edges'

    edges_line = [(i, text) for (i, text) in enumerate(pajek_text) if text.lower() == "*edges\n"][0][0]
    for edgeline in pajek_text[(edges_line+1):]:
        i, j = edgeline.rstrip().split(' ')
        G.add_edge(int(i)-1, int(j)-1)

    print (nx.info(G))
    assert nx.is_bipartite(G)

    #node set level not graph level: print ('Average bipartite clustering:', nx.bipartite.average_clustering(G))
    #node level not graph level: print ('Latapy clustering:', nx.bipartite.latapy_clustering(G)
    print ('Robins-Alexander bipartite clustering: ',  bipartite.robins_alexander_clustering(G))
    # node level not graph level: print ('Square clustering:', nx.square_clustering(G))




    
if __name__ == "__main__":
    main()

