#!/usr/bin/env python
##############################################################################
#
# snowballSampleFromEdgeList.py - snowball sample from edge list
#
# File:    snowballSampleFromEdgeList.py
# Author:  Alex Stivala
# Created: May 2018
#
#
##############################################################################
### TODO update to use EstimNetDirected output format and Pajek input format
"""Do snowball sampling in a (large) network, retaining zone information
 for each sampled node.

 Input file is csv with two columns giving edge list of graph: each
 column is a node id 1..N e.g.:

 1,2
 1,3
 ...etc.

 The graph may be directed or undirected. If directed we do snowball
 sampling on the unidrected version of the graph (i.e. ignore edge
 directions), and the sampled graph is the directed subgraph of the
 original directed graph induced by the nodes thus sampled.

 Output files (sample description file giving names of following files,
 subgraphs as dense matrices, zone files giving zone for each node,
 attirbute files giving attributes for each node) in a directory
 in format used by parallel SPNet.

 WARNING: the output files are overwritten if they exist.

 Reimplementation of the R/igraph version using Python and SNAP
 instead as R is too slow and not scalable.
 
 For SNAP see

 http://snap.stanford.edu/snappy/index.html

 Used version 4.0.0.

 Usage:
 
 python snowballSampleFromEdgelist.py [-d]  adjlist.csv num_samples num_seeds num_waves outputdir

    -d : graph is directed, otherwise undirected
    adjlist.csv is .csv file with edge list as above
    num_samples is number of snowball samples to create
    num_seeds it number of seeds in each sample
    num_Waves is numer of snowball sampling waves
    outputdir is output directory to create output files in

"""

import sys,os,time
import getopt
import random

import snap

from snowballSample import snowball_sample,write_graph_file,write_zone_file,write_subactors_file


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
    sys.stderr.write("usage: " + progname + " [-d] adjlist.csv num_samples num_seeds num_waves outputdir\n")
    sys.stderr.write("  -d: directed graph\n")
    sys.exit(1)


def main():
    """
    See usage message in module header block
    """
    directed = False
    try:
        opts,args = getopt.getopt(sys.argv[1:], "d")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        if opt == '-d':
            directed = True
        else:
            usage(sys.argv[0])

    if len(args) != 5:
        usage(sys.argv[0])

    edgelistFilename = args[0]
    num_samples = int(args[1])
    num_seeds = int(args[2])
    num_waves = int(args[3]) - 1 # -1 for consistency with SPNet
    outputdir = args[4]

    print "directed:", directed
    print "number of samples:", num_samples
    print "number of seeds:", num_seeds
    print "number of waves:", num_waves
    print "output directory:", outputdir
    
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    G = snap.LoadEdgeList(snap.PNGraph if directed else snap.PUNGraph, 
                          edgelistFilename, 0, 1, ',')
    snap.PrintInfo(G)


    # get num_samples * num_seeds distinct random seed nodes (sample without replacement)
    # and convert to list of lists where each list is seed set for one sample
    allseeds = random.sample([node.GetId() for node in G.Nodes()], num_samples * num_seeds)
    seedsets = [allseeds[i:i+num_seeds] for i in range(0, len(allseeds), num_seeds)]

    sampledesc_filename = outputdir + os.path.sep + "sampledesc" + os.path.extsep + "txt"
    sampledesc_f = open(sampledesc_filename, 'w')

    for i in range(num_samples):
        sys.stdout.write( 'generating snowball sample ' + str(i+1) + '... ' )
        start = time.time()
        # have to convert seedset to TIntV for SNAP
        seedsVec = snap.TIntV()
        for nodeid in seedsets[i]:
            seedsVec.Add(nodeid)
        Gsample0 = snowball_sample(G, num_waves, seedsVec)
#        print 'XXX',Gsample0.GetIntAttrDatN(Gsample0.GetRndNId(), "zone")#XXX
        # renumber nodes so they are numbered 0..N-1
        # Actually can't do this as it loses the node attributes (zone)
        # so instead build a dictionary mapping nodeid:zone and use it to get zone
        # in same order as nodes written in adjacency matrix on output
        #Gsample = snap.ConvertGraph(snap.PNEANet, Gsample0, True)
#        print 'YYY',Gsample.GetIntAttrDatN(Gsample.GetRndNId(), "zone")#XXX
        Gsample = Gsample0
        nodelist = list()  # keep this iteration in list so we always use same order in future
        zonedict = dict() # map nodeid : zone
        for node in Gsample.Nodes():
            nodelist.append(node.GetId())
            zonedict[node.GetId()] = Gsample.GetIntAttrDatN(node.GetId(), "zone")
        print time.time() - start, 's'
        
        snap.PrintInfo(Gsample)
        subgraph_filename = outputdir + os.path.sep + "subgraph" + str(i) + os.path.extsep + "txt"
        write_graph_file(subgraph_filename, Gsample, nodelist)
        subzone_filename = outputdir + os.path.sep + "subzone" + str(i) + os.path.extsep + "clu"
        write_zone_file(subzone_filename, Gsample, nodelist, zonedict)
        subactor_filename = outputdir + os.path.sep + "subactor" + str(i) + os.path.extsep + "txt"
        # TODO get actor attributes (currently just writes file with no attrs)
        write_subactors_file(subactor_filename, Gsample, nodelist)
        
        # format of sampledesc file is:
        # N subzone_filename subgraph_filename subactor_filename
        sampledesc_filename = outputdir + os.path.sep + "sampledesc" + os.path.extsep + "txt"
        sampledesc_f.write("%d %s %s %s\n" % (Gsample.GetNodes(), subzone_filename,
                                              subgraph_filename, subactor_filename))

    sampledesc_f.close()

        

    
if __name__ == "__main__":
    main()


