#!/usr/bin/env python
#
# File:    changeStatisticsDirected.py
# Author:  Maksym Byshkin, Alex Stivala
# Created: October 2017
#
#
#
# Functions to compute change statistics.  Each function takes a
# Digraph G (which contains graph itself, along with binary
# attributes and two-path matrices for fast computation) and returns
# the change statistic for adding the edge i,j
#
# These functions are adapted from the original PNet code by Peng Wang:
#
#   Wang P, Robins G, Pattison P. PNet: A program for the simulation and
#   estimation of exponential random graph models. University of
#   Melbourne. 2006.
#

import math
import numpy as np         # used for matrix & vector data types and functions

from Digraph import Digraph

lambda0 = 2.0  #decay factor for alternating statistics (lambda a python keyword)



def changeArc(G, i ,j):
    """
    change statistic for Arc
    """
    return 1

def changeReciprocity(G, i, j):
    """
    change statistic for Reciprocity
    """
    return int(G.isArc(j, i))

def changeAltInStars(G, i, j):
    """
    change statistic for alternating k-in-stars (popularity spread, AinS)
    """
    assert(lambda0 > 1)
    jindegree = G.indegree(j)
    return lambda0 * (1 - (1-1/lambda0)**jindegree)

def changeAltOutStars(G, i, j):
    """
    change statistic for alternating k-out-stars (activity spread, AoutS)
    """
    assert(lambda0 > 1)
    ioutdegree = G.outdegree(i)
    return lambda0 * (1 - (1-1/lambda0)**ioutdegree)

def changeAltKTrianglesT(G, i, j):
    """
    change statistic for alternating k-triangles AT-T (path closure)
    """
    assert(lambda0 > 1)
    delta = 0
    for v in G.outIterator(i):
        assert(G.isArc(i, v))
        if v == i or v == j:
            continue
        if G.isArc(j, v):
            delta += (1-1/lambda0)**G.MixTwoPathMatrix[i, v]
    for v in G.inIterator(i):
        assert(G.isArc(v, i))
        if v == i or v == j:
            continue
        if G.isArc(v, j):
            delta += (1-1/lambda0)**G.MixTwoPathMatrix[v, j]
    delta += lambda0 * (1 - (1-1/lambda0)**G.MixTwoPathMatrix[i, j])
    return delta                 

def changeAltKTrianglesC(G, i, j):
    """
    change statistic for alternating k-triangles AT-C (cyclic closure)
    """
    assert(lambda0 > 1)
    delta = 0
    for v in G.inIterator(i):
        assert(G.isArc(v, i))
        if v == i or v == j:
            continue
        if G.isArc(j, v):
            delta += ( (1-1/lambda0)**G.MixTwoPathMatrix[i, v] +
                       (1-1/lambda0)**G.MixTwoPathMatrix[v, j] )
    delta += lambda0 * (1 - (1-1/lambda0)**G.MixTwoPathMatrix[j, i])
    return delta                 

def changeAltTwoPathsT(G, i, j):
    """
    change statistic for alternating two-paths A2P-T (multiple 2-paths)
    """
    assert(lambda0 > 1)
    delta = 0
    for v in G.outIterator(j):
        assert(G.isArc(j, v))
        if v == i or v == j:
            continue
        delta += (1-1/lambda0)**G.MixTwoPathMatrix[i, v]
    for v in G.inIterator(i):
        assert(G.isArc(v, i))
        if v == i or v == j:
            continue
        delta += (1-1/lambda0)**G.MixTwoPathMatrix[v, j]
    return delta

def changeAltTwoPathsD(G, i, j):
    """
    change statistic for alternating two-paths A2P-D (shared popularity)
    """
    assert(lambda0 > 1)
    delta = 0
    for v in G.outIterator(i):
        assert(G.isArc(i, v))
        if v == i or v == j:
            continue
        delta += (1-1/lambda0)**G.OutTwoPathMatrix[j, v]
    return delta

def changeAltTwoPathsTD(G, i, j):
    """
    change statistic for altnernating two-paths A2P-TD (shared
    popularity + multiple two-paths), adjusting for multiple counting
    """
    return 0.5*(changeAltTwoPathsT(G, i, j) + changeAltTwoPathsD(G, i, j))

def changeSender(G, i, j):
    """
    change statistic for Sender
    """
    return G.binattr[i]

def changeReceiver(G, i, j):
    """
    change statistic for Receiver
    """
    return G.binattr[j]

def changeInteraction(G, i, j):
    """
    change statistic for Interaction
    """
    return G.binattr[i] * G.binattr[j]

def changeMatching(G, i, j):
    """
    change statistic for categorical matching
    """
    return G.catattr[i] == G.catattr[j]

def changeMatchingReciprocity(G, i, j):
    """
    change statistic for categorical matching reciprocity
    """
    return G.catattr[i] == G.catattr[j] and G.isArc(j, i)


def changeMismatching(G, i, j):
    """
    change statistic for categorical mismatching
    """
    return G.catattr[i] != G.catattr[j]

def changeMismatchingReciprocity(G, i, j):
    """
    change statistic for categorical mismatching reciprocity
    """
    return G.catattr[i] != G.catattr[j] and G.isArc(j, i)
