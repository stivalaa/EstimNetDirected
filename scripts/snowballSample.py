##############################################################################
#
# snowballSample.py -  snowball sampling in networks
#
#
# File:    snowballSample.py
# Author:  Alex Stivala
# Created: April 2018
#
##############################################################################
"""
Functions to do snowball sampling in a (large) network using 
SNAP, annotating node with snowball sampling zone.

Reimplementation of the R/igraph version using Python and SNAP
instead as R is too slow and not scalable.

For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.0.0.
"""


import snap


#-----------------------------------------------------------------------------
#
# Functions
#
#-----------------------------------------------------------------------------

def snowball_sample(G, num_waves, seeds):
    """
    Parameters:
      G - SNAP graph or network to sample frpm
      num_waves - number of snowball waves 
      seeds - SNAP vector (TIntV) of seeds (node ids) to start snowball sample 
             from
    
    Return value:
      SNAP network (TNEANet) snowball sampled from G with each node having 
      an integer "zone" attribute for snowball sampling zone 
       (0=seed, 1=first wave, etc.)
      [TNEANet needed to allow zone attribute, not actually using multigraph 
       capability].

    Note directions on directed graph are ignored - can sample in undirected
    or directed graph.
    """
    assert(len(seeds) == len(set(seeds))) # no duplicate node ids
    # It seems like GetSubGraph does not preserve node attributse
    # so instead of adding attributes ot nodes on N, make a Python
    # dictionary mapping node ids to zone and then add them back
    # ass attributes on the subgraph (node ids are preserved so we
    # can do this)
    zonedict = dict() # map nodeid : zone
    N = snap.ConvertGraph(snap.PNEANet, G) # copy graph/network G to network N
    nodes = set(seeds) # will accumulate all nodes (including seeds) here
    for seed in seeds:
        zonedict[seed] = 0  # seed nodes are zone 0
    newNodes = set(nodes)
    for i in range(num_waves):
        wave = i + 1
        #print 'wave',wave
        for node in set(newNodes):
            neighbours = snap.TIntV()
            snap.GetNodesAtHop(G, node, 1, neighbours, False) # neighbours of node
            newNeighbours = set(neighbours) - nodes # neighbours that are not already in nodes
            for node in newNeighbours:
                if not zonedict.has_key(node):
                    zonedict[node] = wave
            newNodes.update(newNeighbours) # newNodes gets set union of itslf and newNeighbours
        nodes.update(newNodes)
    # have to convert nodes set into TIntV for use in SNAP 
    NodeVec = snap.TIntV()
    for node in nodes:
        NodeVec.Add(node)
    sampleN = snap.GetSubGraph(N, NodeVec)
    # now put the zones as attributes on the subgraph nodes (which depends
    # on nodeids being preserved in the subgraph)
    sampleN.AddIntAttrN("zone", -1)  # add zone attribute init to -1
    for (nodeid, zone) in zonedict.iteritems():
        sampleN.AddIntAttrDatN(nodeid, zone, "zone")
    return sampleN



def write_graph_file(filename, G, nodelist, write_header=True):
    """
    write_graph_file() - write edge list in Pajek format
    
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
    with open(filename, 'w') as f:
        if write_header:
            f.write("*vertices " + str(G.GetNodes()) + "\n")
            f.write("*arcs\n")
        for EI in G.Edges():
            f.write("%d %d\n" % (seqdict[EI.GetSrcNId()], seqdict[EI.GetDstNId()]))


def write_zone_file(filename, G, nodelist, zonedict):
    """
    write_zone_file() - write zone file in EstimNetDirected format
    
    The format of the zone file is just the header line "zone"
    and the the zone (staring at 0) of each node one per line.

    Parameters:
      filename -filename to write to (warning: overwritten)
      G - SNAP graph/network object. Must
           have "zone" attribute on nodes 
      nodelist - list of nodeids used to order the nodes in the output
      zonedict - dictionary mapping nodeid to zone
          
    Return value:
      None
    """
    assert(len(zonedict) == len(nodelist))
    assert(len(nodelist) == G.GetNodes())
    with open(filename, 'w') as f:
        #f.write("*vertices " + str(G.GetNodes()) + '\n')
        f.write("zone\n")
        for i in nodelist:
            assert(G.GetIntAttrDatN(i, "zone") == zonedict[i])
            f.write(str(G.GetIntAttrDatN(i, "zone")) + '\n')


