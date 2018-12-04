#!/usr/bin/env python
##############################################################################
#
# convertNBERpatentDataToEstimNetdirectedFormat.py - convert NBER patent data
#
# File:    convertNBERpatentDataToEstimNetdirectedFormat.py
# Author:  Alex Stivala
# Created: December 2018
#
#
##############################################################################

"""Convert NBER patent and citation data to EstimNetDirected format.

Output files (sample description file giving names of following files,
subgraphs as Pajek edge lists (node numbers 1..N_s), zone files giving
zone for each node, attirbute files giving attributes for each node)
in a directory in format used by EstimNetDirected.

Usage:
 
   convertNBERpatentDataToEstimNetdirectedFormat.py data_dir 

   data_dir is directory containing the patent citation data from NBER
                         http://www.nber.org/patents/

 Output files in cwd (WARNING overwritten):
     patent_citations.txt
     patent_binattr.txt
     patent_catattr.txt
     patent_contattr.txt

 WARNING: the output files are overwritten if they exist.

  See

  http://www.nber.org/patents/pat63_99.txt

  for description of variables

References:

Hall, B., Jaffe, A., & Trajtenberg, M. (2001). The NBER patent
citations data file: Lessons, insights and methodological tools. NBER
working paper no. 8498.

Jaffe, A. B., & Trajtenberg, M. (2002). Patents, citations, and
innovations: A window on the knowledge economy. MIT press.


For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.1.0.

NB this uses at least 5 GB memory and tmp directory space

"""

import sys,os,time
import getopt
import random
import math

import snap

from load_nber_patent_data import load_nber_patent_data
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


# def write_attributes_file_binary(filename, G, nodelist, patdata, colnames):
#     """
#     write_attributes_file_binary() - write binary node attribute file 
    
#     The EstimNetDirected format of the binary actor attribute file is 
#     the header line with attribute names and then
#     the attribute value for each on one line per node.  See
#     load_integer_attributes() in digraph.c

#     Parameters:
#         filename -filename to write to (warning: overwritten)
#         G - SNAP graph/network object.
#         nodelist - list of nodeids used to order the nodes in the output
#         patdata - dictionary mapping node ID (int) to list
#                   of attributes (all strings)
#         colnames - dict mapping attribute name to 
#                   index of the patdata list so e.g. we can look
#                   up APPYEAR of patent id 123 with 
#                    patdata[123][colnames['APPYEAR']]
          
#     Return value:
#       None
#     """
#     assert(len(nodelist) == G.GetNodes())
#     assert(len(patdata) >= G.GetNodes())
#     binattrs = ['gender', 'public']
#     # rename gender to male for binary attribute
#     binattr_names = ['male' if x == 'gender' else x for x in binattrs] 
#     with open(filename, 'w') as f:
#         f.write(' '.join(binattr_names) + '\n')
#         for i in nodelist:
#             for attr in binattrs:
#                 val = patdata[i][colnames[attr]]
#                 val = val if val in ['0','1'] else 'NA'
#                 f.write(val)
#                 if attr == binattrs[-1]:
#                     f.write('\n')
#                 else:
#                     f.write(' ' )


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
    catattrs = ['COUNTRY', 'POSTATE', 'ASSIGNEE', 'ASSCODE', 'NCLASS',  #USPTO original variables
                'CAT', 'SUBCAT']   # constructed variables
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
    assert(len(nodelist) == G.GetNodes())
    assert(len(patdata) >= G.GetNodes())
    contattrs = ['APPYEAR', 'GYEAR', 'GDATE', 'CLAIMS']
    contattr_names = contattrs
    with open(filename, 'w') as f:
        f.write(' '.join(contattr_names) + '\n')
        for i in nodelist:
            for attr in contattrs:
                val = patdata[i][colnames[attr]]
                if val != "NA":
                    if val == '':
                        val = "NA"
                    else:
                        val = float(val)
                        if math.isnan(val):
                            val = "NA"
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
    sys.stderr.write("usage: " + progname + " data_dir\n")
    sys.exit(1)


def main():
    """
    See usage message in module header block
    """
    directed = True
    try:
        opts,args = getopt.getopt(sys.argv[1:], "")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        usage(sys.argv[0])

    if len(args) != 1:
        usage(sys.argv[0])

    data_dir = args[0]

    outputdir = '.'

    sys.stdout.write('loading data from ' + data_dir + '...')
    start = time.time()
    (G, patdata, colnames) = load_nber_patent_data(data_dir)
    print time.time() - start, 's'

    snap.PrintInfo(G)


    # We do not add attributes to nodes as SNAP node attribute as
    # these seem to get lost by varoius operations including subgraph
    # that we need to use, so instead maintain them just in the
    # dictionary mapping the original node ids to the attributes -
    # fortunately the original node ids are maintained by
    # GetSubGraph() so we can used these to index the patdata
    # dictoinary in the subgraphs


    # Cannot do this:
    #patdata[:][colnames['COUNTRY']] = convert_to_int_cat(patdata[:][colnames['COUNTRY']]) # like factor in R
    # as get "TypeError: unhashable type" so have to do this instead:
    id_countries = [(k, p[colnames['COUNTRY']]) for (k,p) in patdata.iteritems()]
    id_countries_int = convert_to_int_cat([x[1] for x in id_countries])
    for i in xrange(len(id_countries)):
        patdata[id_countries[i][0]][colnames['COUNTRY']] = id_countries_int[i]
    for attr in ['COUNTRY']:
        sys.stdout.write('There are %d NA for %s\n' % ([p[colnames[attr]] for p in patdata.itervalues()].count('NA'), attr))

    id_states = [(k, p[colnames['POSTATE']]) for (k,p) in patdata.iteritems()]
    id_states_int = convert_to_int_cat([x[1] for x in id_states])
    for i in xrange(len(id_states)):
        patdata[id_states[i][0]][colnames['POSTATE']] = id_states_int[i]
    for attr in ['POSTATE']:
        sys.stdout.write('There are %d NA for %s\n' % ([p[colnames[attr]] for p in patdata.itervalues()].count('NA'), attr))


    # There are 3774768 unique patent identifiers in the citation data but
    # only 2923922 unique patent identifiers in the patent data (patdata).
    # The size of the set intersection of these patent ids is 2755865
    # i.e. there is patent data for 73% of the patents in the citation network.
    # Presumably this is because the patdata (pat63_99.txt) contains all
    # utilit patents in the period 1963 to 1999 but the citation data
    # cit75_99.txt contains all US patent citations for utility patents
    # granted in the period 1975 to 1999, so there are patent ids in here
    # from earlier periods that are cited by patents in that period,
    # for which therefore we don't have the patent data (prior to 1963).
    # So we have to set the data for all patents in network that we have it
    # for, and the rest (27%) to NA.


    citepatent_count = 0
    patentdata_count = 0
    nodelist = list()  # keep this iteration in list so we always use same order in future
    for node in G.Nodes():
        citepatent_count += 1
        patid = node.GetId()
        nodelist.append(patid)
        print citepatent_count, patentdata_count, patid  #XXX
        if not patdata.has_key(patid):
            print 'NA for ', patid #XXX
            patdata[patid] = len(colnames)*["NA"]
        else:
            patentdata_count += 1
    sys.stdout.write("There are %d unique cited/citing patents of which %d (%f%%) have patent data\n" % (citepatent_count, patentdata_count, 100*float(patentdata_count)/citepatent_count))

    graph_filename = outputdir + os.path.sep + "patent_citations" + os.path.extsep + "txt"
    write_graph_file(graph_filename, G, nodelist)
    attributes_binary_filename = outputdir + os.path.sep + "patent_binattr" + os.path.extsep + "txt"
    attributes_categorical_filename = outputdir + os.path.sep + "patent_catattr"  + os.path.extsep + "txt"
    attributes_continuous_filename = outputdir + os.path.sep + "patent_contattr" + os.path.extsep + "txt"

    # write_attributes_file_binary(attributes_binary_filename, G, nodelist, citpatdata, colnames)
    write_attributes_file_categorical(attributes_categorical_filename, G, nodelist, patdata, colnames)
    write_attributes_file_continuous(attributes_continuous_filename, G, nodelist, patdata, colnames)

    nodeid_filename = outputdir + os.path.sep + "subnodeid" + str(i) + os.path.extsep + "txt"
    write_subgraph_nodeids(nodeid_filename, nodelist)

        

    
if __name__ == "__main__":
    main()


