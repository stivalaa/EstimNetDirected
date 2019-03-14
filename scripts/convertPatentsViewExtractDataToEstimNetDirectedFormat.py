#!/usr/bin/env python
##############################################################################
#
# convertPatentsViewExtractDataToEstimNetdirectedFormat.py - 
#          convert PatentsView extracted data to EstimNetDirected format
#
# File:    convertPatentsViewExtractDataToEstimNetdirectedFormat.py
# Author:  Alex Stivala
# Created: March 2019
#
#
##############################################################################

"""
Convert extracts PatentsView patent and citation data 
to EstimNetDirected format.

Output files (sample description file giving names of following files,
subgraphs as Pajek edge lists (node numbers 1..N_s), zone files giving
zone for each node, attirbute files giving attributes for each node)
in a directory in format used by EstimNetDirected.

Usage:
 
   convertPatentsViewExtractDataToEstimNetdirectedFormat.py [-d] data_dir 

   data_dir is directory containing the extracted patent citation data

   -d : only use subgraph of patents that have attribute data



 Output files in cwd (WARNING overwritten):
     patent_citations.txt
     patent_binattr.txt   [currently not used]
     patent_catattr.txt
     patent_contattr.txt
     nodeid.txt

 WARNING: the output files are overwritten if they exist.

For PatentsView see

http://www.patentsview.org/query

For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.1.0.

NB this uses around 7 GB memory and 5 GB tmp directory space

"""

import sys,os,time
import getopt
import math
from datetime import datetime 

import snap

from load_patentsview_extract_data import load_patentsview_extract_data,patch_years
from snowballSample import write_graph_file


#-----------------------------------------------------------------------------
#
# Functions
#
#-----------------------------------------------------------------------------


def convert_to_int_cat(attrs):
    """
    convert_to_int_cat() - convert string categorical attrs to integer

    Like factor() in R, convert categories represented as strings into
    integers.

    Parameters:
       attrs - list of string attributes
    
    Return value:
       list of integer attributes corresponding to the strings
    
    """
    # build dict mapping string to integer for unique strings in attrs list
    fdict = dict([(y,x) for (x,y) in enumerate(set(attrs))])
    print(fdict) # output for possible future reversal (TODO write to file)
    return ['NA' if x == '' else fdict[x] for x in attrs]


def str_to_float(sval):
    """
    str_to_float() - convert string to floating point handling NA

    Parameters:
        sval - string to convert

    Return value:
        floating point value or "NA"
    """
    val = sval
    if val != "NA":
        if val == '':
            val = "NA"
        else:
            val = float(val)
            if math.isnan(val):
                val = "NA"
    return val



def write_attributes_file_categorical(filename, G, nodelist, patdata, colnames):
    """
    write_attributes_file_categorical() - write categorical node attribute file 
    
    The EstimNetDirected format of the categorical actor attribute file is 
    the header line with attribute names and then
    bthe attribute value (integer) for each on one line per node.  See
    load_integer_attributes() in digraph.c

    Parameters:
        filename -filename to write to (warning: overwritten)
        G - SNAP graph/network object.
        nodelist - list of nodeids used to order the nodes in the output
        patdata - dictionary mapping node ID (int) to list
                  of attributes (all strings)
        colnames - dict mapping attribute name to 
                  index of the patdata list so e.g. we can look
                  up APPYEAR of patent id 123 with 
                   patdata[123][colnames['APPYEAR']]
          
    Return value:
      None
    """
    assert(len(nodelist) == G.GetNodes())
    assert(len(patdata) >= G.GetNodes())
    catattrs = ['filing_country', "techcategory_nber", "techsubcategory_nber"]
    catattr_names = catattrs
    with open(filename, 'w') as f:
        f.write(' '.join(catattr_names) + '\n')
        for i in nodelist:
            for attr in catattrs:
                val = patdata[i][colnames[attr]]
                val = val if isinstance(val, int) else (int(val) if val.isdigit() else 'NA')
                f.write(str(val))
                if attr == catattrs[-1]:
                    f.write('\n')
                else:
                    f.write(' ' )


def write_attributes_file_continuous(filename, G, nodelist, patdata, colnames):
    """
    write_attributes_file_continuous() - write continuous node attribute file 
    
    The EstimNetDirected format of the continuous actor attribute file is 
    the header line with attribute names and then
    bthe attribute value (integer) for each on one line per node.  See
    load_integer_attributes() in digraph.c

    Parameters:
        filename -filename to write to (warning: overwritten)
        G - SNAP graph/network object.
        nodelist - list of nodeids used to order the nodes in the output
        patdata - dictionary mapping node ID (int) to list
                  of attributes (all strings)
        colnames - dict mapping attribute name to 
                  index of the patdata list so e.g. we can look
                  up APPYEAR of patent id 123 with 
                   patdata[123][colnames['APPYEAR']]
          
    Return value:
      None
    """

    ## convert dates:
    ##   datetime.strptime('1969-6-39', '%Y-%m-%d').toordinal()
    ##  "proleptic Gregorian ordinal"

    assert(len(nodelist) == G.GetNodes())
    assert(len(patdata) >= G.GetNodes())
    contattrs = ['grantdate','num_claims','filing_date']
    contattr_names = contattrs
    with open(filename, 'w') as f:
        f.write(' '.join(contattr_names) + '\n')
        for i in nodelist:
            for attr in contattrs:
                if attr == 'grantdate' or attr == 'filing_date':
                    if patdata[i][colnames[attr]] != "NA":
                        try:
                            val = datetime.strptime(patdata[i][colnames[attr]],
                                               '%Y-%m-%d').toordinal()
                        except ValueError:
                            val = "NA"
                else:
                    val = str_to_float(patdata[i][colnames[attr]])
                f.write(str(val))
                if attr == contattrs[-1]:
                    f.write('\n')
                else:
                    f.write(' ' )




def write_subgraph_nodeids(filename, nodelist):
    """write_subgraph_nodeids() - write mapping from subgraph sequential ids
                              to original graph node ids

    Writes the original graph node identifiers in file one per line in
    same order as zones and attributes so we can cross-reference the
    subgraph nodes back to the original grpah if necessary.  First
    line is just header "nodeid" than next line is original node id of
    node 1 in subgraph, etc.

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
    sys.stderr.write("usage: " + progname + "  data_dir\n"
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
    (G, patdata, colnames) = load_patentsview_extract_data(data_dir)
    print time.time() - start, 's'

    snap.PrintInfo(G)

    # Remove loops (self-edges).
    # There is actually for some reason one loop (patent id 5489070).
    # G is a PNGraph so multiple edges not allowed in this type anyway.
    snap.DelSelfEdges(G)
    snap.PrintInfo(G)

    # We do not add attributes to nodes as SNAP node attribute as
    # these seem to get lost by varoius operations including subgraph
    # that we need to use, so instead maintain them just in the
    # dictionary mapping the original node ids to the attributes -
    # fortunately the original node ids are maintained by
    # GetSubGraph() so we can used these to index the patdata
    # dictoinary in the subgraphs

    # columns are:
    # patent_id,grantdate,num_claims,filing_country,filing_date,techcategory_nber,techsubcategory_nber
    # e.g.:
    # 10000000,2018-06-19,20,US,2015-03-10,,
    # 3930277,1976-01-06,7,US,1974-08-21,6,69
    # 9999999,2018-06-19,2,US,2015-12-07,,

    id_countries = [(k, p[colnames['filing_country']]) for (k,p) in patdata.iteritems()]
    id_countries_int = convert_to_int_cat([x[1] for x in id_countries])
    for i in xrange(len(id_countries)):
        patdata[id_countries[i][0]][colnames['filing_country']] = id_countries_int[i]
    for attr in ['filing_country']:
        sys.stdout.write('There are %d NA for %s\n' % ([p[colnames[attr]] for p in patdata.itervalues()].count('NA'), attr))


    nodelist = list()  # keep the iteration below in list so we always use same order in future

    if get_subgraph:
        # get subgraph induced by nodes that have patent data 
        nodeVec = snap.TIntV() # nodelist in TIntV format for use in SNAP
        for node in G.Nodes():
            patid = node.GetId()
            if patdata.has_key(patid):
                nodelist.append(patid)
                nodeVec.Add(patid)
        G = snap.GetSubGraph(G, nodeVec)
        print 'Subgraph with only nodes with patent attribute data:'
        snap.PrintInfo(G)
    else:
        # keep all the graph and just put NA for all data attributes
        citepatent_count = 0
        patentdata_count = 0
        for node in G.Nodes():
            citepatent_count += 1
            patid = node.GetId()
            nodelist.append(patid)
            #print citepatent_count, patentdata_count, patid  #XXX
            if not patdata.has_key(patid):
                #print 'NA for ', patid #XXX
                patdata[patid] = len(colnames)*["NA"]
            else:
                patentdata_count += 1
        sys.stdout.write("There are %d unique cited/citing patents of which %d (%f%%) have patent data\n" % (citepatent_count, patentdata_count, 100*float(patentdata_count)/citepatent_count))

    patch_years(data_dir, patdata, colnames) # fix up some wrong date data

    graph_filename = outputdir + os.path.sep + "patent_citations" + os.path.extsep + "txt"
    write_graph_file(graph_filename, G, nodelist)
    attributes_categorical_filename = outputdir + os.path.sep + "patent_catattr"  + os.path.extsep + "txt"
    attributes_continuous_filename = outputdir + os.path.sep + "patent_contattr" + os.path.extsep + "txt"

    write_attributes_file_categorical(attributes_categorical_filename, G, nodelist, patdata, colnames)
    write_attributes_file_continuous(attributes_continuous_filename, G, nodelist, patdata, colnames)

    nodeid_filename = outputdir + os.path.sep + "nodeid" + os.path.extsep + "txt"
    write_subgraph_nodeids(nodeid_filename, nodelist)

        

    
if __name__ == "__main__":
    main()




