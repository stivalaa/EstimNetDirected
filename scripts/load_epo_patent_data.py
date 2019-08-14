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
        # it really means col 1 is citing patent, col 2 is cited patent
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
    # Note separator is '|' not ',' as ',' is used as delimiter of list of classes
    # e.g.
    # EP0000001|1978|DE|LU|FISW GMBH|F25,B27,G06,F28,B23,G11,H04
    # but PatID column 0 used as dict key so skip it
    # and we remove the leading 'EP' on each patent to get (still unique) integer ids
    colnames = csviter.next()[1:] # skip PatID column 0
    # append new column nanmes for data added later
    newcolnames = ['NumClasses','English','Switzerland','Belgium','Sections','NumSections', 'SectionA', 'SectionB','SectionC', 'SectionD', 'SectionE', 'SectionF', 'SectionG', 'SectionH', 'SectionY','CPCsections', 'LanguageCode', 'CountryCode','French','German']
    colnames += newcolnames
    patent_colnames = dict([(name, col) for (col, name) in enumerate(colnames)])
    # have already read header line so rest of iterable csv read is the data
    patentdata = [ (x[0][2:], x[1:] + len(newcolnames)*['NA']) for x in  csviter] 
    patentdict = dict([(int(x[0]), x[1]) for x in patentdata])

    # Now add the new column NumClasses which is just length of
    # the comma-delimited list of classes
    for patid in patentdict.iterkeys():
        patentdict[patid][patent_colnames['NumClasses']] =  len(patentdict[patid][patent_colnames['Classes']].split(','))

    # Build binary attributes for languages and country Switzerland
    for patid in patentdict.iterkeys():
        patentdict[patid][patent_colnames['English']] = ('NA' if patentdict[patid][patent_colnames['Language']] == '' or patentdict[patid][patent_colnames['Language']] == 'XX' else  (1 if patentdict[patid][patent_colnames['Language']] == 'EN' else 0))
        patentdict[patid][patent_colnames['French']] = ('NA' if patentdict[patid][patent_colnames['Language']] == '' or patentdict[patid][patent_colnames['Language']] == 'XX' else  (1 if patentdict[patid][patent_colnames['Language']] == 'FR' else 0))
        patentdict[patid][patent_colnames['German']] = ('NA' if patentdict[patid][patent_colnames['Language']] == '' or patentdict[patid][patent_colnames['Language']] == 'XX' else  (1 if patentdict[patid][patent_colnames['Language']] == 'GR' else 0))
        patentdict[patid][patent_colnames['Switzerland']] = ('NA' if patentdict[patid][patent_colnames['Country']] == '' else (1 if patentdict[patid][patent_colnames['Country']] == 'CH' else 0))
        patentdict[patid][patent_colnames['Belgium']] = ('NA' if patentdict[patid][patent_colnames['Country']] == '' else (1 if patentdict[patid][patent_colnames['Country']] == 'BE' else 0))

    # Add CPC section which is the first character (letter) in the class
    # and also number of sections each patent is in.
    # There are 9 CPC sections:
    # https://www.epo.org/searching-for-patents/helpful-resources/first-time-here/classification/cpc.html
    # A = Human necessities
    # B = Performing operations; transporting 
    # C = Chemistry; metallurgy 
    # D = Textiles; paper
    # E = Fixed constructions
    # F = Mechanical engineering; lighting; heating; weapons; blasting engines or pumps 
    # G = Physics
    # H = Electicity
    # Y = General tagging of new technological developments; general tagging of cross-sectional technologies spanning over several sections of the IPC; technical subjects covered by former USPC cross-reference art collections [XRACs] and digests
    for patid in patentdict.iterkeys():
        # Sections and CPCsections are the same, but Sections will be 
        # converted to int set later and CPCsections will not;
        # we will only list unique CPCsections here (but leave Sections
        # as the R 'factor' like code will do that later)
        patentdict[patid][patent_colnames['Sections']] =  ','.join([tclass[0] for tclass in patentdict[patid][patent_colnames['Classes']].split(',')])
        patentdict[patid][patent_colnames['CPCsections']] =  ','.join(set([tclass[0] for tclass in patentdict[patid][patent_colnames['Classes']].split(',')]))

    # Now add the new column NumSections which is number of unique sections in
    # the comma-delimited list of sections
    # note we use set here to make sure we count unique sections
    for patid in patentdict.iterkeys():
        patentdict[patid][patent_colnames['NumSections']] =  len(set(patentdict[patid][patent_colnames['Sections']].split(',')))

    # Build binary attribute for each CPC Section
    for patid in patentdict.iterkeys():
        for cpcsection in ['A','B','C','D','E','F','G','H','Y']:
            patentdict[patid][patent_colnames['Section'+cpcsection]] = (1 if cpcsection in patentdict[patid][patent_colnames['Sections']].split(',') else 0)
    
    # LanguageCode and CountryCode are jsut Language and Country
    # but kept as original strings for convenience in R not made integers
    # for EstimNetDirected
    for patid in patentdict.iterkeys():
        patentdict[patid][patent_colnames['CountryCode']] = patentdict[patid][patent_colnames['Country']]
        patentdict[patid][patent_colnames['LanguageCode']] = patentdict[patid][patent_colnames['Language']]
        
    
    return (G, patentdict, patent_colnames)
