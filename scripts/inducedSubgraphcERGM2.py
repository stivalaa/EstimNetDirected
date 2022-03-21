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
##############################################################################
"""
Functions to get induced subgraph for citation ERGM (cERGM) 
goodness-of-fit tests in a (large) network using 
SNAP, annotating node with term (time period)

Reimplementation of the R/igraph code in plotEstimNetDirectedcERGMSimFit2.r
using Python and SNAP instead as R is too slow and not scalable: 
the original R/igraph code takes several hours per network on the patent
data (approx 100 000 nodes per time period) for example, making it
unusuable since large numbers of networks are needed (and the cluster has
a 48 hour time limit, and even embarrassing parallelism is not practical either
since in R even 'multicore' means the entire memory is copied for each
core),  while this takes approx 0.1 seconds, and also uses less memory.



For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.1.0.

NB Using Python 2.7 (not Python 3) as could not get SNAP to install on Python 3.
"""

import gzip

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
      G - SNAP network to take subgraph of, must have 'term' int attribute
    
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
    maxterm = max([G.GetIntAttrDatN(i, "term") for i in G.Nodes()])
    maxterm_nodes = [node.GetId() for node in G.Nodes() if G.GetIntAttrDatN(node, "term") == maxterm]
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
    subgraphG = snap.GetSubGraph(G, NodeVec)
    # now put the terms as attributes on the subgraph nodes (which depends
    # on nodeids being preserved in the subgraph)
    subgraphG.AddIntAttrN("term", -1)  # add term attribute init to -1
    for (nodeid, term) in termdict.iteritems():
        subgraphG.AddIntAttrDatN(nodeid, term, "term")
    return subgraphG



def write_term_file_compressed(filename, G, nodelist, termdict):
    """
    write_term_file() - write term file in EstimNetDirected format gzipped
    
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
    with gzip.open(filename, 'w') as f:
        f.write("term\n".encode())
        for i in nodelist:
            assert(G.GetIntAttrDatN(i, "term") == termdict[i])
            f.write((str(G.GetIntAttrDatN(i, "term")) + '\n').encode())



def write_graph_file_compressed(filename, G, nodelist, write_header=True):
    """
    write_graph_file_compressed() - write edge list in Pajek format gzipped
    
    Note that because snap.ConvertGraph() fails to keep node
    attributes so we cannot use it to renumber nodes, we also
    use nodelist to get a sequential node number for each node:
    the index in nodelist of each nodeid is its sequential number,
    so we can these out as sequential node numbers 1..N

    So in order to do this we have to write the output ourselves in
    an iterator, cannot use the snap.SavePajek() function.

    Parameters:
      filename - filename to write to (warning: overwritten)
      G - SNAP graph object
      nodelist - list of nodeids, used to order the nodes in the output
      write_header - if True write Pajek header lines

    Return value:
       None
    """
    assert(len(nodelist) == G.GetNodes())
    assert(len(nodelist) == len(set(nodelist))) # nodeids must be unique
    # build dict mapping nodeid to sequential node number 1..N
    seqdict = {nodeid:(seq+1) for seq, nodeid in enumerate(nodelist)}
    with gzip.open(filename, 'wb') as f:
        if write_header:
            f.write("*vertices " + (str(G.GetNodes()) + "\n").encode())
            f.write("*arcs\n".encode())
        for EI in G.Edges():
            f.write(("%d %d\n" % (seqdict[EI.GetSrcNId()], seqdict[EI.GetDstNId()])).encode())


