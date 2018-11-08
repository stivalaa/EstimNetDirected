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
    G = load_pokec_data('/home/stivala/SNAPestimations/pokec/')

NB this uses at least 0.5 GB memory and tmp directory space

"""

import os,sys
import glob
import tempfile
import gzip
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
    




def load_pokec_data(indirname):
    """Load the pokec data from specified directory

    
    Parameters:
       indirname - path name of directory to load from

    Return value:
       tuple(G, profile) where
        G - SNAP TNGraph object built from the data
        profile - dictionary mapping node ID (int) to list
                  of attributes (all strings)
        profile_colnames - dict mapping attribute name to 
                  index of the profile list so e.g. we can look
                  up AGE of userid 123 with 
                   profile[123][profile_colnames['AGE']]

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
        fin = gzip.open(os.path.join(indirname, infilename), 'rb')
        filename = os.path.join(tmpdir, "soc-pokec-relationships.txt")
        fout = open(filename, 'w')
        fout.write(fin.read())
	fout.close()
        G = snap.LoadEdgeList(snap.PNGraph, filename, 0, 1, '\t')
    finally:
        cleanup_tmpdir(tmpdir)

    # https://snap.stanford.edu/data/soc-pokec-readme.txt
    # but 'user_id' column 0 used as dict key so not included here
    colnames = [         'public', 'completion_percentage',
                'gender', 'region', 'last_login', 'registration',
                'AGE', 'body', 'I_am_working_in_field',
                'spoken_languages', 'hobbies',
                'I_most_enjoy_good_food', 'pets', 'body_type',
                'my_eyesight', 'eye_color', 'hair_color',
                'hair_type', 'completed_level_of_education',
                'favourite_color', 'relation_to_smoking',
                'relation_to_alcohol', 'sign_in_zodiac',
                'on_pokec_i_am_looking_for', 'love_is_for_me',
                'relation_to_casual_sex', 'my_partner_should_be',
                'marital_status', 'children',
                'relation_to_children', 'I_like_movies',
                'I_like_watching_movie', 'I_like_music',
                'I_mostly_like_listening_to_music',
                'the_idea_of_good_evening',
                'I_like_specialties_from_kitchen', 'fun',
                'I_am_going_to_concerts', 'my_active_sports',
                'my_passive_sports', 'profession', 'I_like_books',
                'life_style', 'music', 'cars', 'politics',
                'relationships', 'art_culture',
                'hobbies_interests', 'science_technologies',
                'computers_internet', 'education', 'sport',
                'movies', 'travelling', 'health',
                'companies_brands', 'more']
    profile_colnames = dict([(name, col) for (col, name) in enumerate(colnames)])
    profilepath = os.path.join(indirname, "soc-pokec-profiles.txt.gz")
    profiledata = [ (x[0], x[1:]) for x in csv.reader(gzip.open(profilepath, 'rb'), delimiter='\t') ]
    profiledict = dict([(int(x[0]), x[1]) for x in profiledata])
    assert(G.GetNodes() == len(profiledict))
    return (G, profiledict, profile_colnames)
