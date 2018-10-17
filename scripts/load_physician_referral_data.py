##############################################################################
#
# load_physician_referral_data.py -  load physician referral network in SNAP
#
#
# File:    load_physician_referral_data.py
# Author:  Alex Stivala
# Created: April 2018
#
##############################################################################

"""
Function to load the physician referral data from zip file and
convert to SNAP format.

Physician referral data from Centers for Medicare & Medicaid Services
(CMS.gov): https://questions.cms.gov/faq.php?faqId=7977
(that address was used on 13 Dec 2017 when I downloaded the data, now
seems to be located at:
https://www.cms.gov/Regulations-and-Guidance/Legislation/FOIA/Referral-Data-FAQs.html).
  
as used in the paper:
  
 An, C., O'Malley, A. J., Rockmore, D. N., & Stock, C. D. (2017).
 Analysis of the US patient referral network. Statistics in Medicine. 
 DOI: 10.1002/sim.7565
  
Documentation from 
http://downloads.cms.gov/foia/physician_shared_patient_patterns_technical_requirements.pdf
  
 As in the paper we use the 30 day interval data. Also 2014 as the most
 recent complete year (2015 only has first 7 months):
 http://downloads.cms.gov/foia/physician-shared-patient-patterns-2014-days30.zip

For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.0.0.

E.g. 
    G = load_physician_referral_data('/vlsci/VR0261/stivalaa/Physician_referral_data/physician-shared-patient-patterns-2014-days30.zip')

NB this uses at least 2.5 GB tmp file space and memory

"""

import os,sys
import glob
import tempfile
import zipfile

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
    


def load_physician_referral_data(infilename):
    """
    Load the US physician referral data from specified zipfile

    
    Parameters:
       infilename - path name of zipflie to load from

    Return value:
       SNAP TNGraph object built from the data
    """
    tmpdir = tempfile.mkdtemp()
    try:
        archive = zipfile.ZipFile(infilename, 'r')
        archive.extract('physician-shared-patient-patterns-2014-days30.txt', tmpdir)
        filename = os.path.join(tmpdir, "physician-shared-patient-patterns-2014-days30.txt")
        archive.close()
        context = snap.TTableContext()
        schema = snap.Schema()
        ## schema.Add(snap.TStrTAttrPr("NPI_1", snap.atInt))
        ## schema.Add(snap.TStrTAttrPr("NPI_2", snap.atInt))
        # the above 2 lines worked with SNAP 4.0.0 on VLSCI 
        # but now using SNAP 4.1.0
        # on hpc.ics.usi.ch find that all ids are -1 so graph wrong.
        # Cannot work out why so changed to string not int to try to fix it:
        schema.Add(snap.TStrTAttrPr("NPI_1", snap.atStr))
        schema.Add(snap.TStrTAttrPr("NPI_2", snap.atStr))
        ## schema.Add(snap.TStrTAttrPr("count", snap.atInt))
        ## schema.Add(snap.TStrTAttrPr("unique_bene", snap.atInt))
        ## schema.Add(snap.TStrTAttrPr("same_day_count", snap.atInt))
        # The above 3 lines also worked fine with SNAP 4.0.0 before but
        # now fail on SNAP 4.1.0 (seems to be due to spaces in CSV fields,
        # not inexplicable like first two which have no spaces) but not using
        # them at the moment anyway so easier to just make (unused) strings:
        schema.Add(snap.TStrTAttrPr("count", snap.atStr))
        schema.Add(snap.TStrTAttrPr("unique_bene", snap.atStr))
        schema.Add(snap.TStrTAttrPr("same_day_count", snap.atStr))
        table = snap.TTable.LoadSS(schema, filename, context, ",", snap.TBool(False))
        G = snap.ToGraph(snap.PNGraph, table, "NPI_1", "NPI_2", snap.aaFirst) 
    finally:
        cleanup_tmpdir(tmpdir)

    return G


