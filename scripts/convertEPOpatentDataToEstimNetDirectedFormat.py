#!/usr/bin/env python
##############################################################################
#
# convertEPOpatentDataToEstimNetdirectedFormat.py - convert EPO patent data
#
# File:    convertEPOpatentDataToEstimNetdirectedFormat.py
# Author:  Alex Stivala
# Created: March 2019
#
#
##############################################################################

"""Convert EPO patent and citation data to EstimNetDirected format.

Output files (sample description file giving names of following files,
subgraphs as Pajek edge lists (node numbers 1..N_s), zone files giving
zone for each node, attirbute files giving attributes for each node)
in a directory in format used by EstimNetDirected.

Usage:
 
   convertEPOpatentDataToEstimNetdirectedFormat.py [-d] data_dir 

   data_dir is directory containing the patent citation data from EPO
            extract at USI Informatics
                         
   -d : only use subgraph of patents that have attribute data


 Output files in cwd (WARNING overwritten):
     patent_citations.txt
     patent_binattr.txt
     patent_catattr.txt
     patent_contattr.txt
     patent_setattr.txt
     nodeid.txt
     patent_string_categories.txt [not used by EstimNetDirected but useful to read in R]

 WARNING: the output files are overwritten if they exist.

For SNAP see

http://snap.stanford.edu/snappy/index.html

Used version 4.1.0.

"""

import sys,os,time
import getopt
import random
import math

import snap

from load_epo_patent_data import load_epo_patent_data
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
    return ['NA' if (x == '' or x == 'XX') else fdict[x] for x in attrs]


def convert_to_int_set(attrs):
    """
    convert_to_int_cat_set() - convert string categorical list attrs to 
                               integer sets

    Parametrs:
        attrs - list of list of string attributes

    Return value:
        list of integer sets corresponding to the string lists
    """
    allstrings = [item for sublist in attrs for item in sublist]
    fdict = dict([(y,x) for (x,y) in enumerate(set(allstrings))])
    print(fdict) # output for possible future reversal (TODO write to file)
    return map(lambda x : set([fdict[y] for y in x]), attrs)



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


def write_attributes_file_binary(filename, G, nodelist, patdata, colnames):
    """
    write_attributes_file_binary() - write categorical node attribute file 
    
    The EstimNetDirected format of the binary actor attribute file is 
    the header line with attribute names and then
    bthe attribute value (0 or 1) for each on one line per node.  See
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
    binattrs = ['English','Switzerland','Belgium','SectionA', 'SectionB','SectionC', 'SectionD', 'SectionE', 'SectionF', 'SectionG', 'SectionH', 'SectionY','French','German']
    with open(filename, 'w') as f:
        f.write(' '.join(binattrs) + '\n')
        for i in nodelist:
            for attr in binattrs:
                val = patdata[i][colnames[attr]]
                val = str(val) if val in [0,1] else 'NA'
                f.write(val)
                if attr == binattrs[-1]:
                    f.write('\n')
                else:
                    f.write(' ' )



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
    catattrs = ['Language','Country']
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
    contattrs = ['NumClasses',    # number of technology classes
                 'NumSections',   # number of technology sections (highest level)
                 'Year',          # in data: application year
                 'YearBase1978' ] # constructed here
    contattr_names = contattrs
    with open(filename, 'w') as f:
        f.write(' '.join(contattr_names) + '\n')
        for i in nodelist:
            for attr in contattrs:
                if attr == 'YearBase1978':
                    # year - 1978
                    # for sender or receiver year effects
                    val = str_to_float(patdata[i][colnames["Year"]])
                    if val != "NA":
                        val = val - 1978 if val >= 1978 else "NA"
                else:
                    val = str_to_float(patdata[i][colnames[attr]])
                f.write(str(val))
                if attr == contattrs[-1]:
                    f.write('\n')
                else:
                    f.write(' ' )


def write_attributes_file_set(filename, G, nodelist, patdata, colnames):
    """
    write_attributes_file_set() - write set node attribute file 
    
    The EstimNetDirected format of the set actor attribute file is 
    the header line with attribute names and then
    the comma-delimited list of (integer categories) in the set on
    each line for each attribute.

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
    setattrs = ['Classes', 'Sections']    #  technology classes and sections
    setattr_names = setattrs
    with open(filename, 'w') as f:
        f.write(' '.join(setattr_names) + '\n')
        for i in nodelist:
            for attr in setattrs:
                if attr == 'Classes' or attr == 'Sections':
                    val = patdata[i][colnames[attr]]
                else:
                    assert(False)
                if val is None:
                    f.write('NA')
                elif len(val) == 0:
                    f.write('none')
                else:
                    f.write(','.join([str(x) for x in val]))
                if attr == setattrs[-1]:
                    f.write('\n')
                else:
                    f.write(' ' )




def write_subgraph_nodeids(filename, nodelist):
    """write_subgraph_nodeids() - write mapping from subgraph sequential ids
                              to original graph node ids

    Writes the original graph node identifiers in file one per line in
    same order as zones and attributes so we can cross-reference the
    subgraph nodes back to the original grpah if necessary.  First
    line is just header "nodeid" than next line is original 
    patent identifier (nodeid) of node 1 in subgraph, etc.

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
    sys.stderr.write("usage: " + progname + " [-d] data_dir\n"
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
    (G, patdata, colnames) = load_epo_patent_data(data_dir)
    print time.time() - start, 's'

    snap.PrintInfo(G)

    # Remove loops (self-edges).
    # There is actually for some reason 92 nodes with self-loops
    # e.g. EP0021443
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


    # convert categorical attribute values to integers like factor in R
    for cat_colname in ['Language','Country']:
        catvalues = [(k, p[colnames[cat_colname]]) for (k,p) in patdata.iteritems()]
        catvalues_int = convert_to_int_cat([x[1] for x in catvalues])
        for i in xrange(len(catvalues)):
            patdata[catvalues[i][0]][colnames[cat_colname]] = catvalues_int[i]
        sys.stdout.write('There are %d NA for %s\n' % ([p[colnames[cat_colname]] for p in patdata.itervalues()].count('NA'), cat_colname))


    # convert categorical set attribute values to integers like factor in R
    for set_colname in ['Classes','Sections']:
        setvalues = [(k, p[colnames[set_colname]]) for (k,p) in patdata.iteritems()]
        setvalues_int = convert_to_int_set([x[1].split(',') for x in setvalues])
        for i in xrange(len(setvalues)):
            patdata[setvalues[i][0]][colnames[set_colname]] = setvalues_int[i]
        sys.stdout.write('There are %d NA for %s\n' % ([p[colnames[set_colname]] for p in patdata.itervalues()].count('NA'), set_colname))

    nodelist = list()  # keep the iteration below in list so we always use same order in future

    if get_subgraph:
        # get subgraph induced by nodes that have patent data in the
        # pat63_99.txt file
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


    graph_filename = outputdir + os.path.sep + "patent_citations" + os.path.extsep + "txt"
    write_graph_file(graph_filename, G, nodelist)
    attributes_binary_filename = outputdir + os.path.sep + "patent_binattr"  + os.path.extsep + "txt"
    attributes_categorical_filename = outputdir + os.path.sep + "patent_catattr"  + os.path.extsep + "txt"
    attributes_continuous_filename = outputdir + os.path.sep + "patent_contattr" + os.path.extsep + "txt"
    attributes_set_filename = outputdir + os.path.sep + "patent_setattr" + os.path.extsep + "txt"

    write_attributes_file_binary(attributes_binary_filename, G, nodelist, patdata, colnames)
    write_attributes_file_categorical(attributes_categorical_filename, G, nodelist, patdata, colnames)
    write_attributes_file_continuous(attributes_continuous_filename, G, nodelist, patdata, colnames)
    write_attributes_file_set(attributes_set_filename, G, nodelist, patdata, colnames)

    nodeid_filename = outputdir + os.path.sep + "nodeid" + os.path.extsep + "txt"
    write_subgraph_nodeids(nodeid_filename, nodelist)

    # write patent sections as original letters before converting to int
    # This cannot be used by EstimNetDirected but is useful to read in R
    # and factor there so that the original names are preserved
    sections_filename = outputdir + os.path.sep + "patent_string_categories" + os.path.extsep + "txt"
    attrnames = ['CPCsections','LanguageCode','CountryCode']
    with open(sections_filename, 'w') as f:
        f.write(' '.join(attrnames) + '\n')
        for i in nodelist:
            for attrname in attrnames:
                val = patdata[i][colnames[attrname]]
                val = 'NA' if (val == 'NA' or val == 'XX') else val
                f.write(val)
                if attrname == attrnames[-1]:
                    f.write('\n')
                else:
                    f.write(' ' )

    
if __name__ == "__main__":
    main()

