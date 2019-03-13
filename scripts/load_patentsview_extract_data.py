##############################################################################
#
# load_patentsview_extract_data.py -  load extracted PatentsView data
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
    patentsview_patch.zip: 
             patent_utility_patchednberyear.csv      - corrected years for some patents
             own_cat_subcat_patentview.csv           - recoding of category and subcategory

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
    # e.g.:
    # 10000000,2018-06-19,20,US,2015-03-10,,
    # 3930277,1976-01-06,7,US,1974-08-21,6,69
    # 9999999,2018-06-19,2,US,2015-12-07,,
    # but patent_id column 0 used as dict key so skip it
    colnames = csviter.next()[1:] # skip patent_id column 0
    patent_colnames = dict([(name, col) for (col, name) in enumerate(colnames)])
    # have already read header line so rest of iterable csv read is the data
    # remove the patents that have an id that is not an int (starts with a 
    # letter) as these are not utility patents (design patents etc.)
    patentdata = [ (x[0], x[1:]) for x in  csviter]
    patentdict = dict([(int(x[0]), x[1]) for x in patentdata if x[0].isdigit()])
    return (G, patentdict, patent_colnames)



def patch_years(indirname, patentdict, patent_colnames):
    """Patch the years on some entries with bad dates

    
    Parameters:
        indirname - path name of directory to load from
        patentdict - dictionary mapping patent ID (int) to list
                  of attributes (all strings) [IN/OUT]
        patent_colnames - dict mapping attribute name to 
                  index of the patent list so e.g. we can look
                  up filing_date of userid 123 with 
                   patent[123][patent_colnames['filing_date']]

    Return value:
        None. Updates patentdict.
    """
    patchzip = os.path.join(indirname, "patentsview_patch.zip")
    zf = zipfile.ZipFile(patchzip)
    csviter = csv.reader(zf.open("patent_utility_patchednberyear.csv"))
    #  skip header line patent_id,appyear,gyear
    # Then data in each row is e.g.:
    # 3933359,1974,1976
    csviter.next() # skip header
    for (patid, appyear, gyear) in csviter:
        patid = int(patid)
        if patentdict.has_key(patid):
            print 'XXX patch_years',patid,appyear,gyear
            patentdict[patid][patent_colnames['filing_date']] = appyear + patentdict[patid][patent_colnames['filing_date']][4:]
            patentdict[patid][patent_colnames['grantdate']] = gyear + patentdict[patid][patent_colnames['grantdate']][4:]
