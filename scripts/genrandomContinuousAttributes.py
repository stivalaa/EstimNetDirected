#!/usr/bin/env python
##############################################################################
#
# genrandomContinuousAttributes.py - generate random continuous attributes
#
# Generates noramlly distributed continuous attributes with specified
# mean and standard deviation (default mu=0 sigma=1)
#
# File:    genrandomContinuousAttributes.py
# Author:  Alex Stivala
# Created: July 2014
#
##############################################################################

"""
Generates noramlly distributed continuous attributes with specified
mean and standard deviation (default mu=0 sigma=1)

Usage:
    genrandomContinuousAttributes.py n [mu sigma]
    
    n is the size of the sample
    mu is mean
    sigma is standard deviation
   
Output is to stdout in format to be read by PNet for simulation etc.
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
    sys.stderr.write("usage: " + progname + " n [mu sigma]\n")
    sys.stderr.write("   n is the size of the sample\n");
    sys.stderr.write("   mu is the mean(default 0)\n");
    sys.stderr.write("   sigma is the standard deviation(default 1)\n");
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

    if len(args) != 1 and len(args) != 3:
        usage(sys.argv[0])

    mu =0
    sigma = 1
    n = int(args[0])
    if len(args) > 1:
      mu = float(args[1])
      sigma = float(args[2])

    sys.stdout.write("continuousAttribute\n")
    for i in xrange(n):
      sys.stdout.write("%f\n" %(random.normalvariate(mu, sigma)))

    
if __name__ == "__main__":
    main()


