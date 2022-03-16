##############################################################################
#
# inducedSubgraphcERGM2.py -  get subgraphs for cERGM Gof plots
#
#
# File:    inducedSbugGraphcERGM2.py
# Author:  Alex Stivala
# Created: March 2022
#
#
# Status: NOT TESTED YET
#
##############################################################################
"""
Functions to do snowball sampling in a (large) network using 
SNAP, annotating node with term (time period)

Reimplementation of the R/igraph code in plotEstimNetDirectedcERGMSimFit2.r
 using Python and SNAP instead as R is too slow and not scalable.

For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.1.0.

NB Using Python 2.7 (not Python 3) as could not get SNAP to install on Python 3.
"""


import snap


#-----------------------------------------------------------------------------
#
# Functions
#
#-----------------------------------------------------------------------------

def cERGM2_subgraph(G):
    """
    cERGM2_subgraph() - get subgraph for cERGM goodness-of-fit

    Returns the subgraph induced by the union
    of the nodes in the last term (time period), and all nodes in
    earlier terms that receive arcs from them. 
  

    Parameters:
      G - SNAP graph or network to take subgraph of
    
    Return value:
      SNAP network (TNEANet) subgraph of G with each node having 
      an integer "term" attribute for time period (starts at 0)
       (0=seed, 1=first wave, etc.)
      [TNEANet needed to allow term attribute, not actually using multigraph 
       capability].

    """
    # It seems like GetSubGraph does not preserve node attributse
    # so instead of adding attributes ot nodes on N, make a Python
    # dictionary mapping node ids to term and then add them back
    # ass attributes on the subgraph (node ids are preserved so we
    # can do this)
    termdict = dict() # map nodeid : term
    N = snap.ConvertGraph(snap.PNEANet, G) # copy graph/network G to network N
    maxterm = max([G.GetIntAttrDatN(i, "term") for i in G.Nodes()])
    maxterm_nodes = [i for i in G.nodes() if G.GetIntAttdDatN(i, "term") == maxterm]
    nodes = set(maxterm_nodes) # will accumulate all nodes here
    for i in maxterm_nodes:
        termdict[i] = maxterm
    newNodes = set(nodes)
    for node in set(newNodes):
        neighbours = snap.TIntV()
        snap.GetNodesAtHop(G, node, 1, neighbours, True) # out-neighbours of node
        newNeighbours = set(neighbours) - nodes # neighbours that are not already in nodes
        for node in newNeighbours:
            if not termdict.has_key(node):
                termdict[node] = G.GetIntAttrDatN(node, "term")
        newNodes.update(newNeighbours) # newNodes gets set union of itself and newNeighbours
    nodes.update(newNodes)
    # have to convert nodes set into TIntV for use in SNAP 
    NodeVec = snap.TIntV()
    for node in nodes:
        NodeVec.Add(node)
    subgraphN = snap.GetSubGraph(N, NodeVec)
    # now put the terms as attributes on the subgraph nodes (which depends
    # on nodeids being preserved in the subgraph)
    subgraphN.AddIntAttrN("term", -1)  # add term attribute init to -1
    for (nodeid, term) in termdict.iteritems():
        subgraphN.AddIntAttrDatN(nodeid, term, "term")
    return subgraphN


def write_term_file(filename, G, nodelist, termdict):
    """
    write_term_file() - write term file in EstimNetDirected format
    
    The format of the term file is just the header line "term"
    and the the term (staring at 0) of each node one per line.

    Parameters:
      filename -filename to write to (warning: overwritten)
      G - SNAP graph/network object. Must
           have "term" attribute on nodes 
      nodelist - list of nodeids used to order the nodes in the output
      termdict - dictionary mapping nodeid to term
          
    Return value:
      None
    """
    assert(len(termdict) == len(nodelist))
    assert(len(nodelist) == G.GetNodes())
    with open(filename, 'w') as f:
        f.write("term\n")
        for i in nodelist:
            assert(G.GetIntAttrDatN(i, "term") == termdict[i])
            f.write(str(G.GetIntAttrDatN(i, "term")) + '\n')


