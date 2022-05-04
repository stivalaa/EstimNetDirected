#!/usr/bin/env python
################################################################################
#
# File:    bipartitematrix.py
# Author:  Alex Stivala
# Created: July 2012
#
################################################################################
"""
Functions for handling biadjacency matrices and adjacency matrices
"""

import sys

 
#-------------------------------------------------------------------------------
#
# Function definitions
#
#-------------------------------------------------------------------------------

def read_bipartite_matrix(fh):
    """
    Parse a rectangular n x m 0-1 matrix for a bipartite network (the
    biadjacency matrix) from (open for read) filehandle fh and return
    it is a graph in the python dictionary format representing the
    adjacency lists as dictionary of dictionaries, where for any edge
    v-w G[v][w] is 1 if there is an edge or if not. Note here v is
    0..n-1 and w is 0..m-1 and this is implicitly undirected.


    Parameters:
       fh - open for read filehandle to read from

    Return value:
       tuple (m, G) where m is number of columns and
       G is bipartite dictionary of dictionaries as described above
       (note n is len(g))
       or None on error

    """
    m = None
    n = 0
    G = dict()
    for line in fh:
        G[n] = dict()
        sline = line.split()
        if m == None:
            m = len(sline)
        else:
            if len(sline) != m:
                sys.stderr.write("Error: expecting %d columns but got %d "
                                 " on row %d\n" % (m, len(sline), n+1))
                return None
        for j in xrange(m):
            if sline[j] == "1":
                G[n][j] = 1
            elif sline[j] != "0":
                sys.stderr.write("Error: row %d column %d is not 0 or 1\n" %
                                 (n+1, j+1))
                return None
        n += 1
    return (m,G)


def bipartite_to_adjmatrix(m, bpgraph):
    """
    Input is a bipartite graph as returned by read_bipartite_matrix()
    i.e. each bpgraph[v][w] is 1 iff there is an edge v-w, where v in
    0..m-1 is one set of nodes and w in 0..n-1 is the other disjoint
    set of nodes.  Convert this to the general adjacency list graph
    data structure where there are now m+n nodes and the output
    G[v][w] is 1 iff there is an edge between v and w. Note that now
    the nodes are labelled 0..m+n (we hare relablled them in one set,
    not labelled separately for the two disjoint sets that exist in a
    bipartite graph)

    Parameters:
       m - number of columns in bpgraph
       bpgraph - bipartite graph dict as above

    Return value:
       adjlist dict of dicts for m+n square matrix as above
    """
    G = dict()
    n = len(bpgraph)
    for i in xrange(m+n):
        G[i] = dict()
    for i in xrange(n):
        for j in bpgraph[i].iterkeys():
            G[i][n+j] = 1
            G[n+j][i] = 1
    return G
    


def bipartite_to_json(m, bpgraph, comment=None, attributes=None):
    """
    Input is a bipartite graph as returned by read_bipartite_matrix()
    i.e. each bpgraph[v][w] is 1 iff there is an edge v-w, where v in
    0..m-1 is one set of nodes and w in 0..n-1 is the other disjoint
    set of nodes.  Convert this to the JSON representation of a graph
    used by NetworkX JSON (see http://networkx.lanl.gov) and in the
    BPNet REST web service (see Network.java)
    
    Since a matrix representation has no node identifiers, except implicitly
    the row/column indices, the JSON representatino output will label the 
    nodes 1..m+n to ensure unique node identifiers, using the 'bipartite'
    label in the JSON output to mark each node as to its bipartite class (0/1).

    Parameters:
       m - number of columns in bpgraph
       bpgraph - bipartite graph dict as above
       comment - if not None, string to put in a network 'graph' attribute
                 'comment' tuple as value. Must not contain unescaped quote (")
       attributes - if not None, a list, one per node (rows then columns),
                    each element of which is adictionary of keys and values 
                    (each a string whih must have no double quotes) each of
                    which is added as an attribute to the node

    Return value:
       string containing JSON representation of the bipartite graph
    """
    # easier (but less efficient) to convert to general adjacency matrix first
    adjmatrix = bipartite_to_adjmatrix(m, bpgraph)
    n = len(bpgraph)
    assert(len(adjmatrix) == m+n)
    if attributes != None:
        assert(len(attributes) == m+n)
    json_str = '{"graph": [["description", "bipartite"]'
    if comment != None:
        json_str += ',["comment", "' + comment + '"]'
    json_str += '],'
    json_str += '"nodes": ['
    for i in xrange(n+m):
        if i < n:
            bpclass = '0'
        else:
            bpclass = '1'
        json_str += '{"bipartite": ' + bpclass + ', "id": ' + str(i+1)
        if attributes != None:
            json_str += ','
            count = 0
            for (key,value) in attributes[i].iteritems():
                json_str += '"' + key + '":' + value
                if count < len(attributes[i])-1:
                    json_str += ','
                count += 1
        json_str +=  '}'
        if i < m+n-1:
            json_str += ', '
    json_str += '], "adjacency":'
    json_str += bipartite_to_adjlist_json(m, bpgraph, adjmatrix)
    json_str += '}'
    return json_str



def bipartite_to_adjlist_json(m, bpgraph, adjmatrix=None):
    """
    Input is a bipartite graph as returned by read_bipartite_matrix()
    i.e. each bpgraph[v][w] is 1 iff there is an edge v-w, where v in
    0..m-1 is one set of nodes and w in 0..n-1 is the other disjoint
    set of nodes.  Convert this to the JSON representation of
    the adjacency (same as a bipartite_to_json but ONLY the adjacnecy list).
    
    Since a matrix representation has no node identifiers, except implicitly
    the row/column indices, the JSON representatino output will label the 
    nodes 1..m+n to ensure unique node identifiers, using the 'bipartite'
    label in the JSON output to mark each node as to its bipartite class (0/1).

    Parameters:
       m - number of columns in bpgraph
       bpgraph - bipartite graph dict as above
       adjmatrix - precomputed adjacnecy matrix from bipartite_to_adjmatrix
                  for bpgraph, or None to compute here

    Return value:
       string containing JSON representation of the bipartite graph
    """
    if adjmatrix == None:
        adjmatrix = bipartite_to_adjmatrix(m, bpgraph)
    n = len(bpgraph)
    assert(len(adjmatrix) == m+n)
    json_str = '['
    for i in xrange(len(adjmatrix)):
        json_str += '['
        was_prev_edge = False
        for j in xrange(len(adjmatrix)):
            if adjmatrix[i].has_key(j):
                if was_prev_edge:
                    json_str += ','
                was_prev_edge = True
                json_str += '{"id": ' + str(j+1) + '}'
        json_str += ']'
        if i < len(adjmatrix)-1:
            json_str += ','
    json_str += ']'
    return json_str


def adjlist_to_bipartite(network):
    """
    Convert adjacency list (as per NetworkX json adjacency data format)
    to bipartite matrix in dictionary form.

    Parameters:
        network - bipartite graph in NetorkX adjacency data format.
                  Dictionary with key "nodes" value list of "id":id 
                  and "bipartite":<0|1> key:value
                  pairs and key "adjacency" list of lists where outer list
                  corresponds to nodes list, and each inner list is list
                  of "id":id key:value pairs giving the id of nodes
                  the node for that list position is adjacent to.

        Return value:
          tuple (m, bpgraph) where m is number of columns
          dict bpgraph[v][w] with value 1 if there is a v-w edge
          and v and w are now integers with v in 0..n-1 is biparite=0 nodes
          and w in 0..m-1 is bipartite=1 nodes
          
    """
    n = len([node for node in network["nodes"] if node["bipartite"] == 0])
    m = len([node for node in network["nodes"] if node["bipartite"] == 1])
    assert (m + n == len(network["nodes"]))
    assert (len(network["adjacency"]) == n or len(network["adjacency"]) == m+n)

    # make dict id:i mapping ids to 0-based enumeratino for the two clases
    idMap0 = dict(zip([node["id"] for node in network["nodes"] if node["bipartite"] == 0], xrange(n)))
    idMap1 = dict(zip([node["id"] for node in network["nodes"] if node["bipartite"] == 1], xrange(m)))
                  
    
    bpgraph = dict()
    for i in xrange(n):
        bpgraph[i] = dict()

    listi = 0
    for adjlist in network["adjacency"]:
        snode = network["nodes"][listi]
        if snode["bipartite"] == 0:
            i = idMap0[snode["id"]]
        else:
            i = idMap1[snode["id"]]
        for edge in adjlist:
            nodeid = edge["id"]
            if idMap0.has_key(nodeid):
                bpclass = 0
                j = idMap0[nodeid]
            elif idMap1.has_key(nodeid):
                bpclass = 1
                j = idMap1[nodeid]
            else:
                raise KeyError("no node with id " + nodeid)
            if bpclass == snode["bipartite"]:
                raise ValueError("network is not bipartite")
            if bpclass == 0:
                bpgraph[j][i] = 1
            else:
                bpgraph[i][j] = 1
        listi += 1
            

    return (m, bpgraph)

def write_bipartite_matrix(fh, m, bpgraph):
    """
    Write a rectangular n x m 0-1 matrix for a bipartite network (the
    biadjacency matrix) to (open for write) filehandle fh.

    

    Parameters:
       fh       - open for write filehandle to write to
       m       - number of coumns in bipartite matrix
       bpgraph - graph in the python dictionary format representing
                 the adjacency lists as dictionary of dictionaries,
                 where for any edge v-w G[v][w] is 1 if there is an
                 edge or if not. Note here v is 0..n-1 and w is 0..m-1
                 and this is implicitly undirected.

    Return value:
       None

    """
    for i in xrange(len(bpgraph)):
        for j in xrange(m):
            if bpgraph[i].has_key(j) and bpgraph[i][j] == 1:
                fh.write("1")
            else:
                fh.write("0")
            if j < m-1:
                fh.write(" ")
        fh.write("\n")


# def jsongraph_to_bipartite(graph):
#     """
#     Convert graph in JSONgraph format (see JSONgraph.java)
#     to bipartite matrix in dictionary form.

#     Parameters:
#         graph -   bipartite graph in JSONgraph data format.
#                   Dictionary with key "nodes" value list of dictionaries with
#                   value 
#                       "name":id 
#                   and "attributes": "bipartite":<0|1>
#                   key:value pairs and 
#                   key "edge" list of dictionary with "src":id and "dst":id
#                   giving id of node the edge connects.

#         Return value:
#           tuple (m, bpgraph) where m is number of columns
#           dict bpgraph[v][w] with value 1 if there is a v-w edge
#           and v and w are now integers with v in 0..n-1 is biparite=0 nodes
#           and w in 0..m-1 is bipartite=1 nodes
          
#     """
#     n = len([node for node in graph["nodes"] if node["attributes"]["bipartite"] == 0])
#     m = len([node for node in graph["nodes"] if node["attributes"]["bipartite"] == 1])
#     assert (m + n == len(network["nodes"]))
#     assert (len(network["adjacency"]) == n or len(network["adjacency"]) == m+n)

#     # make dict id:i mapping ids to 0-based enumeratino for the two clases
#     idMap0 = dict(zip([node["name"] for node in graph["nodes"] if node["attributes"]["bipartite"] == 0], xrange(n)))
#     idMap1 = dict(zip([node["name"] for node in graph["nodes"] if node["attributes"]["bipartite"] == 1], xrange(m)))
                  
    
#     bpgraph = dict()
#     for i in xrange(n):
#         bpgraph[i] = dict()

#     for edge in graph["edges"]:
#         if idMap0.has_key(edge["src"]):
#             i = idMap0[edge["src"]]
#             bpclass_i = 0
#         elif idMap1.has_key(edge["src"]):
#             i = idMap1[edge["src"]]
#             bpclass_i = 1
#         else:
#             raise KeyError("no node with id " + edge["src"])
        
#         if idMap0.has_key(edge["dst"]):
#             j = idMap0[edge["dst"]]
#             bpclass_j = 0
#         elif idMap1.has_key(edge["dst"]):
#             j = idMap1[edge["dst"]]
#             bpclass_j = 1
#         else:
#             raise KeyError("no node with id " + edge["src"])
        
#         if bpclass_i == bpclass_j:
#             raise ValueError("graph is not bipartite")

#         if bpclass_i == 0:
#             bpgraph[i][j] = 1
#         else:
#             bpgraph[j][i] = 1

#     return (m, bpgraph)


def read_continuous_bipartite_matrix(fh):
    """
    Parse a rectangular n x m 0-1 matrix for a bipartite network from
    (open for read) filehandle fh and return it is a graph in the
    python dictionary format representing the adjacency lists as
    dictionary of dictionaries, where for any edge v-w G[v][w] is the
    weight of th edge between v and w Note here v is 0..n-1 and w is
    0..m-1 and this is implicitly undirected.


    Parameters:
       fh - open for read filehandle to read from

    Return value:
       tuple (m, G) where m is number of columns and
       G is bipartite dictionary of dictionaries as described above
       (note n is len(g))
       or None on error

    """
    m = None
    n = 0
    G = dict()
    for line in fh:
        G[n] = dict()
        sline = line.split()
        if m == None:
            m = len(sline)
        else:
            if len(sline) != m:
                sys.stderr.write("Error: expecting %d columns but got %d "
                                 " on row %d\n" % (m, len(sline), n+1))
                return None
        for j in xrange(m):
            G[n][j] = float(sline[j])
        n += 1
    return (m,G)

def continuous_bipartite_to_json(m, G):
    """
    Paramters:
      m - number of columns in bipartite matrix
      G - biparite dictionary of dictionaries as returned by 
          read_continuous_bipartite_matrix()
    
    Return value:
      JSON string contianing the continous values for edge between each
      pair (of different bipartite classes), for use in dyadicCovariates
      BPNet JSON input
    
    """
    n = len(list(G.iteritems())) - m
    jsonstr = '"covariates": ['
    for (i, d) in G.iteritems():
        for (j, v) in d.iteritems():
            jsonstr += '{"src":"' + str(i+1) + '", "dst":"' + str(m+n+j+1) + '", "continuousValue":' + str(v) + '},'
    return jsonstr[:-1] + ']\n' # replace trailing comma with ] then newline

def continuous_bipartite_to_jsongraph(m, G, covname):
    """
    Paramters:
      m - number of columns in bipartite matrix
      G - biparite dictionary of dictionaries as returned by 
          read_continuous_bipartite_matrix()
     covname - name of the dyadic covariate the matrix has values for
    
    Return value: 
      JSON string contianing JSONgraph format continous
      values for edge between each pair (of different bipartite
      classes), for use in dyadicCovariatesJSONgraph BPNet JSON input
    
    """
    n = len(list(G.iteritems())) - m
    jsonstr = '"edges": ['
    for (i, d) in G.iteritems():
        for (j, v) in d.iteritems():
            jsonstr += '{"src":"' + str(i+1) + '", "dst":"' + str(m+n+j+1) + '", "attributes":{"' + covname + '":' + str(v) + '}},'
    return jsonstr[:-1] + ']\n' # replace trailing comma with ] then newline
