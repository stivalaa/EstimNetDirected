##############################################################################
#
# load_pokec_data.py -  load Pokec social network data in SNAP
#
#
# File:    load_pokec_data.py
# Author:  Alex Stivala
# Created: November 2018
#
##############################################################################

"""Function to load the Pokec social network data from zip file and
convert to SNAP format.


The Pokec social network data is from SNAP

https://snap.stanford.edu/data/soc-Pokec.html

See documentation in
https://snap.stanford.edu/data/soc-pokec-readme.txt and source
citation:

L. Takac, M. Zabovsky. Data Analysis in Public Social Networks,
International Scientific Conference & International Workshop
Present Day Trends of Innovations, May 2012 Lomza, Poland.
https://snap.stanford.edu/data/soc-pokec.pdf

Reference for SNAP collection of data sets:

@misc{snapnets,
 author       = {Jure Leskovec and Andrej Krevl},
 title        = {{SNAP Datasets}: {Stanford} Large Network Dataset Collection},
 howpublished = {\url{http://snap.stanford.edu/data}},
 month        = jun,
 year         = 2014
}

Input files (in specified directory):
   soc-pokec-profiles.txt.gz
   soc-pokec-relationships.txt.gz

For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.1.0.

E.g. 
    G = load_snap_data('/home/stivala/SNAPestimations/pokec/')

NB this uses at least 2.5 GB memory and tmp directory space

"""

import os,sys
import glob
import tempfile
import gzip

import snap


#-----------------------------------------------------------------------------
#
# Functions
#
#-----------------------------------------------------------------------------

def cleanup_tmpdir(tmpdir):
    """
    Remove a temporary directory and its contents
    Parameters:
       tmpdir - temporary directory to remove
    Return value: None
    """
    try:
        for filename in glob.glob(os.path.join(tmpdir, "*")):
            os.remove(filename)
        os.rmdir(tmpdir)
    except OSError, inst:
        sys.stderr.write('WARNING: could not remove temp files'
                         ' in ' + tmpdir + '\n' + str(inst) + '\n')
    




def load_pokec_data(indirname):
    """Load the pokec data from specified directory

    
    Parameters:
       indirname - path name of directory to load from

    Return value:
       SNAP TNGraph object built from the data

    Note that in SNAP, node IDs are unique integers and do not have to
    be 0..N-1. However EstimNetDirected requires the node ids in the
    Pajek files for its input are numbered 1..N. Fortunately the SNAP
    Pokec network data has nodes numbered 1..N already so we can just
    directly use those as the node ids and not have to do any
    renumbering etc.
    """
    infilename = "soc-pokec-relationships.txt.gz"
    tmpdir = tempfile.mkdtemp()
    try:
        fin = gzip.open(infilename, 'rb')
        filename = os.path.join(tmpdir, "soc-pokec-relationships.txt")
        fout = open(filename, 'w')
        fout.write(fin.read())
        G = snap.LoadEdgeList(snap.TNGraph, filename, 0, 1, '\t')
    finally:
        cleanup_tmpdir(tmpdir)

    return G
