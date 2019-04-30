#!/usr/bin/env python
##############################################################################
#
# File:    convertGersteinYeastDataToEstimNetDirectdFormat.py 
# Author:  Alex Stivala
# Created: April 2019
#
##############################################################################

"""Load the Gerstein lab yeast regulatory network data in adjacency list
format and convert to EstimNetDirected arclist (Pajek) format using igraph.

The output is to files

yeast_arclist.txt
yeast_binattr.txt
yeast_idmap.txt

The binattr file contains 'self' attibute to mark self-regulating genes.

The yeast_idmap.txt has two columns, first is original name, second is 
sequential id used in the pajek format arclist.

WARNING files are overwritten if they exist.


Downloaded data from

http://archive.gersteinlab.org/proj/Hierarchy_Rewiring/PNAS_hier/YRN.txt

Citation for data

H. Yu and M. Gerstein, "Genomic analysis of the hierarchical structure
of regulatory networks." Proc. Natl. Acad. Sci. USA 103(40),
14724-14731 (2006)

"""

import os,sys

import igraph



def load_from_adjacency_lists(fname):
    """
    load_from_adjacency_lists() - load igraph object from adjacency lists file

    Format of file is each line is tab-delimited, first element is source
    vertex name, subseqeuent are destination vertex names directly
    connected to that source

    Parameters:
        fname = name of file to read from

    Return value:
        igraph Graph object built from adjacnecy lists in fname
    """
    adjlists = [l.split('\t') for l in open(fname).read().splitlines()]
    allgenes = list(set([x for lst in adjlists for x in lst])) # unique ids
    g = igraph.Graph(directed=True)
    g.add_vertices(allgenes)
    for l in adjlists:
        g.add_edges([(l[0], v) for v in l[1:]])

    return(g)




#-----------------------------------------------------------------------------
#
# main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname + "\n")
    sys.exit(1)


def main():
    """
    See usage message in module header block
    """
    if (len(sys.argv) != 1):
        usage(sys.argv[0])

    g = load_from_adjacency_lists("YRN.txt")

    print(g.summary())

    # set binary node attribute 'self' to True if node has self-loop
    g.vs['self'] = [v in v.neighbors() for v in g.vs] 

    # remove self loops (and multiple edges, but there are none of these)
    g = g.simplify()

    print(g.summary())

    g.write_pajek('yeast_arclist.txt')

    f = open('yeast_binattr.txt', 'w')
    f.write('self\n')
    f.write('\n'.join([('1' if x else '0') for x in g.vs['self']]))
    f.write('\n')
    f.close()

    f = open('yeast_idmap.txt', 'w')
    f.write('\n'.join(y +'\t'+ str(x) for (x,y) in enumerate(g.vs['name'])))
    f.write('\n')
    f.close()

if __name__ == "__main__":
    main()
