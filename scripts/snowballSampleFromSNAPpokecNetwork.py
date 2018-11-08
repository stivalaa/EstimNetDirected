#!/usr/bin/env python
##############################################################################
#
# snowballSampleFromSNAPpokecNetwork.py - snowball sample Pokec network
#
# File:    snowballSampleFromSNAPpokecNetwork.py
# Author:  Alex Stivala
# Created: November 2018
#
#
##############################################################################

"""Do snowball sampling in a the Pokec online social network (directed)

Output files (sample description file giving names of following files,
subgraphs as Pajek edge lists (node numbers 1..N_s), zone files giving
zone for each node, attirbute files giving attributes for each node)
in a directory in format used by EstimNetDirected.

Usage:
 
   snowballSampleFromSNAPpokecNetwork.py data_dir num_samples num_seeds num_waves outputdir

   data_dir is directory containing the Pokec data from SNAP 
               (https://snap.stanford.edu/data/soc-Pokec.html)
   num_samples is number of snowball samples to create
   num_seeds it number of seeds in each sample
   num_Waves is numer of snowball sampling waves
   outputdir is output directory to create output files in


 WARNING: the output files are overwritten if they exist.

 For SNAP see

 http://snap.stanford.edu/snappy/index.html

 Used version 4.1.0.

"""

import sys,os,time
import getopt
import random

import snap

from load_pokec_data import load_pokec_data
from snowballSample import snowball_sample,write_graph_file,write_zone_file


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
    sys.stderr.write("usage: " + progname + " data_dir num_samples num_seeds num_waves outputdir\n")
    sys.exit(1)


def main():
    """
    See usage message in module header block
    """
    directed = True
    try:
        opts,args = getopt.getopt(sys.argv[1:], "")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        usage(sys.argv[0])

    if len(args) != 5:
        usage(sys.argv[0])

    data_dir = args[0]
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

    sys.stdout.write('loading data from ' + data_dir + '...')
    start = time.time()
    G = load_pokec_data(data_dir)
    print time.time() - start, 's'

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
        Gsample = snowball_sample(G, num_waves, seedsVec)
        nodelist = list()  # keep this iteration in list so we always use same order in future
        zonedict = dict() # map nodeid : zone
        for node in Gsample.Nodes():
            nodelist.append(node.GetId())
            zonedict[node.GetId()] = Gsample.GetIntAttrDatN(node.GetId(), "zone")
        print time.time() - start, 's'
        
        snap.PrintInfo(Gsample)
        subgraph_filename = outputdir + os.path.sep + "subgraph" + str(i) + os.path.extsep + "txt"
        write_graph_file(subgraph_filename, Gsample, nodelist)
        subzone_filename = outputdir + os.path.sep + "subzone" + str(i) + os.path.extsep + "txt"
        write_zone_file(subzone_filename, Gsample, nodelist, zonedict)
        subactor_filename = outputdir + os.path.sep + "subactor" + str(i) + os.path.extsep + "txt"
        # TODO get actor attributes 
        #write_subactors_file(subactor_filename, Gsample, nodelist)
        
        # format of sampledesc file is:
        # N subzone_filename subgraph_filename subactor_filename
        sampledesc_filename = outputdir + os.path.sep + "sampledesc" + os.path.extsep + "txt"
        sampledesc_f.write("%d %s %s %s\n" % (Gsample.GetNodes(), subzone_filename,
                                              subgraph_filename, subactor_filename))

    sampledesc_f.close()

        

    
if __name__ == "__main__":
    main()


