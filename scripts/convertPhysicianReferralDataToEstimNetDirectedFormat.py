#!/usr/bin/env python
##############################################################################
#
# Convert physican referral data to EstimNetDirected format
#
# File:    convertPhysicianReferralDataToEstimNetdirectedFormat.py
# Author:  Alex Stivala
# Created: April 2020
#
#
##############################################################################

"""Convert physician referral network to EstimNetDirected format.

Output files is physician_referral_arclist.txt in cwd.

Usage:
 
   convertPhysicianReferralDataToEstimNetdirectedFormat.py data_dir 

   data_dir is directory containing the physician referral data from https://questions.cms.gov/faq.php?faqId=7977
            (see load_physician_referral_data.R)
                         
 Output files in cwd (WARNING overwritten):
     physician_referral_arclist.txt
     nodeid.txt

 WARNING: the output files are overwritten if they exist.

For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.1.0.

"""

import sys,os,time
import getopt

import snap

from load_physician_referral_data import load_physician_referral_data
from snowballSample import write_graph_file


#-----------------------------------------------------------------------------
#
# Functions
#
#-----------------------------------------------------------------------------


def write_subgraph_nodeids(filename, nodelist):
    """write_subgraph_nodeids() - write mapping from subgraph sequential ids
                              to original graph node ids

    Writes the original graph node identifiers in file one per line in
    same order as zones and attributes so we can cross-reference the
    subgraph nodes back to the original grpah if necessary.  First
    line is just header "nodeid" than next line is original 
    patent identifier (nodeid) of node 1 in subgraph, etc.

    Paramters:
        filename - filename to write to (warning: overwritten)
        nodelist - list of nodeids used to order the nodes in the output

    Return value:
        None.
    """
    with open(filename, 'w') as f:
        f.write('nodeid\n')
        for i in nodelist:
            f.write(str(i) + '\n')


#-----------------------------------------------------------------------------
#
# main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname + " [-d] data_dir\n"
                     "-d : get subgraph with attribute data nodes only\n")
    sys.exit(1)


def main():
    """
    See usage message in module header block
    """
    get_subgraph = False # if True discard nodes without attribute data
    try:
        opts,args = getopt.getopt(sys.argv[1:], "d")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        if opt == "-d":
            get_subgraph = True
        else:
            usage(sys.argv[0])

    if len(args) != 1:
        usage(sys.argv[0])

    data_dir = args[0]

    outputdir = '.'

    sys.stdout.write('loading data from ' + data_dir + '...')
    start = time.time()
    datazipfile = data_dir + os.path.sep + 'physician-shared-patient-patterns-2014-days30.zip'
    G = load_physician_referral_data(datazipfile)
    print time.time() - start, 's'

    snap.PrintInfo(G)


    # Remove loops (self-edges).
    # G is a PNGraph so multiple edges not allowed in this type anyway.
    snap.DelSelfEdges(G)
    snap.PrintInfo(G)

    # specify ordered nodelist to map sequential ids to original ids consistent
    nodelist = [node.GetId() for node in G.Nodes()]

    graph_filename = outputdir + os.path.sep + "physician_referall_arclist" + os.path.extsep + "txt"
    nodeid_filename = outputdir + os.path.sep + "nodeid" + os.path.extsep + "txt"
    write_graph_file(graph_filename, G, nodelist)
    write_subgraph_nodeids(nodeid_filename, nodelist)

     
if __name__ == "__main__":
    main()

