##############################################################################
#
# load_epo_patent_data.py -  load EPO patent citation data
#
#
# File:    load_epo_patent_data.py
# Author:  Alex Stivala
# Created: March 2019
#
##############################################################################

"""Function to load the EPO patent citation data extract from USI Informatics

For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.1.0.

E.g. 
    (G, patdata, colnames) = load_epo_patent_data('/home/stivala/EPO/')

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
    




def load_epo_patent_data(indirname):
    """Load the EPO patent citation data from specified directory

    
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
    infilename = "Archive.zip"
    tmpdir = tempfile.mkdtemp()
    try:
        zf = zipfile.ZipFile(os.path.join(indirname, infilename))
        filename = os.path.join(tmpdir, "valid-ep-citations.csv")
        csviter = csv.reader(zf.open("valid-ep-citations.csv"))
        fout = open(filename, 'w')
        # skip header line CiteeID,CitedID
        (a,b) = csviter.next()
        assert a == "CiteeID"
        assert b == "CitedID"
        # and we remove the leading 'EP' on each patent to get (still unique) integer ids
        for (id1, id2) in csviter:
            fout.write(id1[2:] + ',' + id2[2:] + '\n')
        fout.close()
        G = snap.LoadEdgeList(snap.PNGraph, filename, 0, 1, ',')
    finally:
        cleanup_tmpdir(tmpdir)

    patentpath = os.path.join(indirname, "Archive.zip")
    zf = zipfile.ZipFile(patentpath)
    csviter = csv.reader(zf.open("patents_ep_selfcontained.csv"), delimiter='|')
    #  get header line 
    #    PatID|Year|Language|Country|Applicant|Classes
    # Note separator is ';' not ',' as ',' is used as delimiter of list of classes
    # e.g.
    # EP0000001|1978|DE|LU|FISW GMBH|F25,B27,G06,F28,B23,G11,H04
    # but PatID column 0 used as dict key so skip it
    # and we remove the leading 'EP' on each patent to get (still unique) integer ids
    colnames = csviter.next()[1:] # skip PatID column 0
    # append new column nanmes for data added later
    newcolnames = ['PrimaryClass']
    colnames += newcolnames
    patent_colnames = dict([(name, col) for (col, name) in enumerate(colnames)])
    # have already read header line so rest of iterable csv read is the data
    patentdata = [ (x[0][2:], x[1:] + len(newcolnames)*['NA']) for x in  csviter] 
    patentdict = dict([(int(x[0]), x[1]) for x in patentdata])

    # Now add the new column PrimaryClass which is just the first class
    # in the comma-delimited list of classes
    for patid in patentdict.iterkeys():
        patentdict[patid][patent_colnames['PrimaryClass']] =  patentdict[patid][patent_colnames['Classes']].split(',')[0]
    return (G, patentdict, patent_colnames)
