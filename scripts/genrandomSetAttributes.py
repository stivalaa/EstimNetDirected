#!/usr/bin/env python
##############################################################################
#
# genrandomSetAttributes.py - generate random set attributes
#
# Generates random set attributes in format for EstimNetDirected setAttrFile
#
# File:    genrandomSetAttributes.py
# Author:  Alex Stivala
# Created: July 2019
#
##############################################################################

"""
Generates random set attributes in format for EstimNetDirected setAttrFile

Usage:
    genrandomSetAttributes.py n
    
    n is the size of the sample
    
Output is to stdout in format to be read by EstimNetDirected in setAttrFile
"""

import sys
import getopt
import random


#-----------------------------------------------------------------------------
#
# Functions
#
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#
# main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname + " n\n")
    sys.stderr.write("   n is the size of the sample\n");
    sys.exit(1)

def main():
    """
    See usage message in module header block
    """
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        usage(sys.argv[0])

    if len(args) != 1:
        usage(sys.argv[0])


    n = int(args[0])
    if len(args) > 1:
      mu = float(args[1])
      sigma = float(args[2])

    sys.stdout.write("setAttribute1 setAttribute2\n")
    for i in xrange(n):
		setsizelist = [10, 20]
		for i in xrange(len(setsizelist)):
			if i > 0:
				sys.stdout.write(' ')
			randset = random.sample(set(range(setsizelist[i])), random.randint(0, setsizelist[i]))
			if len(randset) == 0:
				sys.stdout.write('none')
			else:
				sys.stdout.write(','.join([str(x) for x in randset]))
		sys.stdout.write("\n")


    
if __name__ == "__main__":
    main()


