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
    assignee_patentsview.zip                       - data on assignees
    inventor_patentsview.zip                       - data on inventorss
    patent_utility_patentsview.zip                 - data on granted US utility patents
    uspatentcitation_examinercites_patentsview.zip - examiner citations
    uspatentcitation_utility_patentsview.zip       - US utility patent citations
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
    # append columns owncat and ownsubcat for data to be added later
    newcolnames = ['owncat','ownsubcat']
    # append inventor data columns names for data to be added later
    newcolnames += ['inventor_id','inventor_city','inventor_state','inventor_country','inventor_lat','inventor_long']
    colnames += newcolnames
    patent_colnames = dict([(name, col) for (col, name) in enumerate(colnames)])
    # have already read header line so rest of iterable csv read is the data
    # remove the patents that have an id that is not an int (starts with a 
    # letter) as these are not utility patents (design patents etc.)
    patentdata = [ (x[0], x[1:] + len(newcolnames)*['NA']) for x in  csviter] # append NAs for new columns to have values added later
    patentdict = dict([(int(x[0]), x[1]) for x in patentdata if x[0].isdigit()])

    # get the recoded technology category and subcategory data
    # from patentsview_patch.zip: 
    #    own_cat_subcat_patentview.csv
    patentpatchpath = os.path.join(indirname, "patentsview_patch.zip")
    zf = zipfile.ZipFile(patentpatchpath)
    csviter = csv.reader(zf.open("own_cat_subcat_patentview.csv"))
    # get header line
    #patent_id,mainclass_id_current,subclass_id_current,sequence,mainclass_id,subclass_id,owncat,ownsubcat
    # e.g.
    #3930271,2,2/161.4,0,2,2/161A,6,63
    # (but no documentation for most of these - just using owncat,ownsubcat)
    owncolnames = csviter.next()
    assert owncolnames[0] == 'patent_id'
    assert owncolnames[6] == 'owncat'
    assert owncolnames[7] == 'ownsubcat'
    for (patent_id,mainclass_id_current,subclass_id_current,sequence,mainclass_id,subclass_id,owncat,ownsubcat) in csviter:
        patid = int(patent_id)
        if patentdict.has_key(patid):
            patentdict[patid][patent_colnames['owncat']] = owncat
            patentdict[patid][patent_colnames['ownsubcat']] = ownsubcat
    
    # Get the inventor data:
    #Archive:  inventor_patentsview.zip
    #  inflating: inventor_patentsview.csv
    #patent_id,inventor_id,inventor_sequence,location_id,city,state,country,latlong
    #10000000,10000000-1,0,mwfp6gfqjew8,Manhattan Beach,CA,US,33.8847|-118.41
    inventorpath = os.path.join(indirname, "inventor_patentsview.zip")
    zf = zipfile.ZipFile(inventorpath)
    csviter = csv.reader(zf.open("inventor_patentsview.csv"))
    inventorcolnames = csviter.next()
    assert inventorcolnames[0] == 'patent_id'
    assert inventorcolnames[1] == 'inventor_id'
    assert inventorcolnames[4] == 'city'
    assert inventorcolnames[5] == 'state'
    assert inventorcolnames[6] == 'country'
    assert inventorcolnames[7] == 'latlong'
    for (patent_id,inventor_id,inventor_sequence,location_id,city,state,country,latlong) in csviter:
        if not patent_id[0].isdigit():
            continue  # only want utility patents, which id are all digits
        patid = int(patent_id)
        if patentdict.has_key(patid):
            # only use the inventor id with -1 on the end meaning first inventor
            if inventor_id[-2:] == '-1':
                inentor_id = inventor_id[:-2] # and remove the -1 off the end
                patentdict[patid][patent_colnames['inventor_id']] = inventor_id
                patentdict[patid][patent_colnames['inventor_city']] = city
                patentdict[patid][patent_colnames['inventor_state']] = state
                patentdict[patid][patent_colnames['inventor_country']] = country
                try:
                    (lati,longi) = latlong.split('|')
                except ValueError:
                    (lati,longi) = ('NA','NA')
                patentdict[patid][patent_colnames['inventor_lat']] = lati
                patentdict[patid][patent_colnames['inventor_long']] = longi
            

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
            print 'XXX patch_years',patid,appyear,gyear,
            patentdict[patid][patent_colnames['filing_date']] = appyear + patentdict[patid][patent_colnames['filing_date']][4:]
            patentdict[patid][patent_colnames['grantdate']] = gyear + patentdict[patid][patent_colnames['grantdate']][4:]
            print patentdict[patid][patent_colnames['filing_date']], #XXX
            print patentdict[patid][patent_colnames['grantdate']] #XXX



