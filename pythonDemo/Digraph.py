#!/usr/bin/env python
#
# File:    Digraph.py
# Author:  Maksym Byshkin, Alex Stivala
# Created: October 2017
#
#
#
# Defines the directed graph structure Digraph with arc list graph
# representations and fast lookup matrices.
#
# These functions are adapted from the original PNet code by Peng Wang:
#
#   Wang P, Robins G, Pattison P. PNet: A program for the simulation and
#   estimation of exponential random graph models. University of
#   Melbourne. 2006.
#

import numpy as np         # used for matrix & vector data types and functions
import time


class Digraph:
    """
    The network is represented as a dictionary of dictionaries.
    Nodes are indexed by integers 0..n-1. The outermost dictionary has the node
    v as a key, and dictionary as value. Then this dictionary has the neighbours
    of v as the keys, and here simply 1 as value (there are no edge weights).
    This is the structure suggested by David Eppstein (UC Irvine) in
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/117228
    and as noted there it supprts easy modification by edge addition and
    removal, which is required by Algorithm EE.
    So G[i] is a dictionary with k entries, where k is the degree of node i.
    and G[i][j] exists (and has value 1) exactly when j is a neighbour of i.
    And simple operations are:
      outdegree of node i:                   len(G[i])
      does arc i->j exist?:                  j in G[i]                 
      dict where keys are out-neighbours of i:  G[i]
      iterator over out-neighbourhood of i   G[i].iterkeys()
      in-place insert edge i->j:             G[i][j] = 1
      in-place remove edge i->j:             G[i].pop(j) 
      copy of G (not modified when G is):    deepcopy(G)
      number of arcs:              sum([len(v.keys()) for v in G.itervalues()])
    To make these operations simple (no special cases for nodes with no
    neighbours),  the graph is initialized so that G[i] = dict() for
    all 0 <= i < n: dict(zip(range(n), [dict() for i in range(n)]))

    Node attributes are stored in a separate list, which is simply
    indexed by node id i.e. a simple list of the attribute values in
    node id order.
    """

    def __init__(self, pajek_edgelist_filename, binattr_filename=None,
                 catattr_filename=None):
        """
        Construct digraph from Pajek format network and binary attributes.

        Parameters:
            pajek_edgelist_filename - edge list in Pajek format
            binattr_filename  - binary attributes
                                Default None: no binary attributes loaded
            catattr_filename - categorical attributes
                                Default None: no categorical attributes loaded
        """
        self.G = None  # dict of dicts as described above
        self.Grev = None # version with all arcs reversed to get in-neighbours
        self.binattr = None # binary attributes: list by node (int not boolean)
        self.catattr = None # categorical attributes: list by node

        print 'Reading graph and building matrices...',
        start = time.time()
        f =  open(pajek_edgelist_filename)
        l = f.readline() # first line must be e.g. "*vertices 500"
        n = int(l.split()[1])

        # empty graph n nodes        
        self.G = dict(zip(range(n), [dict() for i in range(n)]))
        self.Grev = dict(zip(range(n), [dict() for i in range(n)]))
        self.InTwoPathMatrix =  np.zeros((n, n))
        self.OutTwoPathMatrix = np.zeros((n, n))
        self.MixTwoPathMatrix = np.zeros((n, n))

        while l.rstrip().lower() != "*arcs":
            l = f.readline()
        lsplit = f.readline().split()
        while len(lsplit) == 2:
            (i, j) = map(int, lsplit)
            assert(i >= 1 and i <= n and j >= 1 and j <= n)
            self.insertArc(i-1, j-1)    # input is 1-based but we are 0-based
            lsplit = f.readline().split()

        if binattr_filename is not None:
            self.binattr = map(int, open(binattr_filename).read().split()[1:]) 
            assert(len(self.binattr) == n)

        if catattr_filename is not None:
            self.catattr = map(int, open(catattr_filename).read().split()[1:]) 
            assert(len(self.catattr) == n)

        print time.time() - start, 's'


    def numNodes(self):
        """
        Return number of nodes in digraph
        """
        return len(self.G)
    
    def numArcs(self):
        """
        Return number of arcs in digraph
        """
        return sum([len(v.keys()) for v in self.G.itervalues()])
    
    def density(self):
        """
        Return the digraph density 
        """
        edges = self.numArcs()
        nodes = self.numNodes()
        return float(edges) / float(nodes*(nodes-1))

    def outdegree(self, i):
        """
        Return Out-degree of node i
        """
        return len(self.G[i])

    def indegree(self, i):
        """
        Return In-degree of node i
        """
        return len(self.Grev[i])
        

    def isArc(self, i, j):
        """
        Return True iff arc i -> j in digraph
        """
        return j in self.G[i]

    def outIterator(self, i):
        """
        Return iterator over out-neighbours of i
        """
        return self.G[i].iterkeys()

    def inIterator(self, i):
        """
        Return iterator over in-neighbours of i
        """
        return self.Grev[i].iterkeys()

    def insertArc(self, i, j):
        """
        Insert arc i -> j in place
        """
        self.G[i][j] = 1
        self.Grev[j][i] = 1
        self.updateTwoPathMatrices(i, j, isAdd=True)

    def removeArc(self, i, j):
        """
        Delete arc i -> j in place
        """
        self.G[i].pop(j)
        self.Grev[j].pop(i)
        self.updateTwoPathMatrices(i, j, isAdd=False)

    def updateTwoPathMatrices(self, i, j, isAdd):
        """
        Update the 2-paths matrices used for fast computation of change
        statistics with arc i->j and isAdd for true for add,
        false for delete.
        """
        incval = 1 if isAdd else -1
        n = self.numNodes()
        for v in xrange(n):
            if v == i or v == j:
                continue
            if self.isArc(i, v):
                self.OutTwoPathMatrix[v,j] += incval
                self.OutTwoPathMatrix[j,v] += incval
            if self.isArc(v, j):
                self.InTwoPathMatrix[v,i] += incval
                self.InTwoPathMatrix[i,v] += incval
            self.MixTwoPathMatrix[v,j] += int(self.isArc(v, i))*incval
            self.MixTwoPathMatrix[i,v] += int(self.isArc(j, v))*incval



