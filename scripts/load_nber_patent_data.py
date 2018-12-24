##############################################################################
#
# load_nber_patent_data.py -  load NBER patent citation data
#
#
# File:    load_nber_patent_data.py
# Author:  Alex Stivala
# Created: December 2018
#
##############################################################################

"""Function to load the NBER patent citation data downloaded from

http://www.nber.org/patents/

References:

Hall, B., Jaffe, A., & Trajtenberg, M. (2001). The NBER patent
citations data file: Lessons, insights and methodological tools. NBER
working paper no. 8498.

Jaffe, A. B., & Trajtenberg, M. (2002). Patents, citations, and
innovations: A window on the knowledge economy. MIT press.

Input files (in specified directory):
   acite75_99.zip
   apat63_99.zip

For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.1.0.

E.g. 
    (G, patdata, colnames) = load_nber_patent_data('/home/stivala/patentCitations/')

NB this uses at least 5 GB memory and tmp directory space

"""

import os,sys
import glob
import tempfile
import zipfile
import csv

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
    




def load_nber_patent_data(indirname):
    """Load the NBER patent citation data from specified directory

    
    Parameters:
       indirname - path name of directory to load from

    Return value:
       tuple(G, patentdict, patent_colnames) where
        G - SNAP TNGraph object built from the data
        patentdict - dictionary mapping patent ID (int) to list
                  of attributes (all strings)
        patent_colnames - dict mapping attribute name to 
                  index of the patent list so e.g. we can look
                  up APPYEAR of userid 123 with 
                   patent[123][patent_colnames['APPYEAR']]

    Note that in SNAP, node IDs are unique integers and do not have to
    be 0..N-1. So The patent ids can be used for these identifiers.
    However EstimNetDirected requires the node ids in the Pajek files
    for its input are numbered 1..N, so we will have to do renumbering
    for the output file (EstimNetDirected input file).
    """
    infilename = "acite75_99.zip"
    tmpdir = tempfile.mkdtemp()
    try:
        zf = zipfile.ZipFile(os.path.join(indirname, infilename))
        filename = os.path.join(tmpdir, "cite75_99.txt")
        fout = open(filename, 'w')
        # skip header line "CITING","CITED"
        fout.write('\n'.join(zf.read("cite75_99.txt").splitlines()[1:]))
	fout.close()
        G = snap.LoadEdgeList(snap.PNGraph, filename, 0, 1, ',')
    finally:
        cleanup_tmpdir(tmpdir)

    # http://www.nber.org/patents/pat63_99.txt
    patentpath = os.path.join(indirname, "apat63_99.zip")
    zf = zipfile.ZipFile(patentpath)
    csviter = csv.reader(zf.open("apat63_99.txt"))
    #  get header line ['PATENT', 'GYEAR', 'GDATE', 'APPYEAR', 'COUNTRY', 'POSTATE', 'ASSIGNEE', 'ASSCODE', 'CLAIMS', 'NCLASS', 'CAT', 'SUBCAT', 'CMADE', 'CRECEIVE', 'RATIOCIT', 'GENERAL', 'ORIGINAL', 'FWDAPLAG', 'BCKGTLAG', 'SELFCTUB', 'SELFCTLB', 'SECDUPBD', 'SECDLWBD']
    # but PATENT column 0 used as dict key so skip it
    colnames = csviter.next()[1:] # skip PATENT column 0
    # add column for binary attribute 1 when there is data about the patent 
    # (because it is in the 1963 - 1999 period in pat63_99.txt) else 0 
    # This field will for patents without data will be set to 0
    # in convertNBERpatentDataToEstimNetDirectedFormat.py
    # main when matching up the patent attributes to citations.
    colnames.append('HASDATA')
    patent_colnames = dict([(name, col) for (col, name) in enumerate(colnames)])
    # have already read header line so rest of iterable csv read is the data
    patentdata = [ (x[0], x[1:] + [1] ) for x in  csviter] #append 1 for HASDATA
    patentdict = dict([(int(x[0]), x[1]) for x in patentdata])
    return (G, patentdict, patent_colnames)
