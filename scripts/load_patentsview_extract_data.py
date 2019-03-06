##############################################################################
#
# load_patentsview_extract_data.py -  load NBER patent citation data
#
#
# File:    load_patentsview_extract_data.py
# Author:  Alex Stivala
# Created: March 2019
#
##############################################################################

"""Function to load the patent PatentsView data query extract

Input files (in specified directory):
    assignee_patentsview.zip                         - data on assignees
    inventor_patentsview.zip                         - data on inventorss
    patent_utility_patentsview.zip                   - data on granted US utility patents
    uspatentcitation_examinercites_patentsview.zip   - examiner citations
    uspatentcitation_utility_patentsview.zip         - US utility patent citations

For PatentsView see

http://www.patentsview.org/query

For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.1.0.

E.g. 
    (G, patdata, colnames) = load_patentsview_extract_data('/home/stivala/PatentViewAlfonsExtract')


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
    




def load_patentsview_extract_data(indirname):
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
                  up filing_date of userid 123 with 
                   patent[123][patent_colnames['filing_date']]

    Note that in SNAP, node IDs are unique integers and do not have to
    be 0..N-1. So The patent ids can be used for these identifiers.
    However EstimNetDirected requires the node ids in the Pajek files
    for its input are numbered 1..N, so we will have to do renumbering
    for the output file (EstimNetDirected input file).
    """
    infilename = "uspatentcitation_utility_patentsview.zip"
    tmpdir = tempfile.mkdtemp()
    try:
        zf = zipfile.ZipFile(os.path.join(indirname, infilename))
        filename = os.path.join(tmpdir, "uspatentcitation_utility_patentsview.csv")
        fout = open(filename, 'w')
        # skip header line "citing","cited"
        fout.write('\n'.join(zf.read("uspatentcitation_utility_patentsview.csv").splitlines()[1:]))
	fout.close()
        G = snap.LoadEdgeList(snap.PNGraph, filename, 0, 1, ',')
    finally:
        cleanup_tmpdir(tmpdir)

    patentpath = os.path.join(indirname, "patent_utility_patentsview.zip")
    zf = zipfile.ZipFile(patentpath)
    csviter = csv.reader(zf.open("patent_utility_patentsview.csv"))
    #  get header line patent_id,grantdate,num_claims,filing_country,filing_date,techcategory_nber,techsubcategory_nber
    # but patent_id column 0 used as dict key so skip it
    colnames = csviter.next()[1:] # skip patent_id column 0
    patent_colnames = dict([(name, col) for (col, name) in enumerate(colnames)])
    # have already read header line so rest of iterable csv read is the data
    # remove the patents that have an id that is not an int (starts with a 
    # letter) as these are not utility patents (design patents etc.)
    patentdata = [ (x[0], x[1:]) for x in  csviter]
    patentdict = dict([(int(x[0]), x[1]) for x in patentdata if x[0].isdigit()])
    return (G, patentdict, patent_colnames)
