#!/usr/bin/env python
#
# File:    EstimNetDirectedSimpleDemo.py
# Author:  Maksym Byshkin, Alex Stivala
# Created: August 2017
#
#
# Simple demonstration implementation of the Equilibrium Expectation algorithm
# for estimation Exponential Random Graph Model (ERGM) parameters from
#
#    Byshkin M, Stivala A, Mira A, Robins G, Lomi A 2018 "Fast
#    maximum likelihood estimation via equilibrium expectation for
#    large network data". Scientific Reports 8:11509
#    doi:10.1038/s41598-018-29725-8
#

from itertools import chain
import time
import os
import random
import math
import numpy as np         # used for matrix & vector data types and functions

from Digraph import Digraph
from changeStatisticsDirected import *



#
# Constants
#

sampler_m = 1000   # number of proposals per iteration

THETA_PREFIX = 'theta_values_' # prefix for theta output filename
DZA_PREFIX = 'dzA_values_'     # prefix for dzA output filename


def BasicSampler(G, changestats_func_list, theta, performMove):
    """
    BasicSampler - sample from ERGM distribution with basic sampler,
                   returning estimate of E(Delta_z(x_obs))

    Parameters:
       G                   - Digraph object for graph to estimate
       changestats_func_list  - list of change statistics funcions
       theta               - numpy vector of theta (parameter) values
       performMove         - if True, actually do the MC move,
                             updating the graph and twoPath matrix
                             (otherwise these are not modified)

    Returns:
        acceptance_rate     - sampler acceptance rate
        addChangeStats      - numpy vector of change stats for add moves
        delChangeStats      - numpy vector of change stats for delete moves

    Note G is updated in place if performMove is True
    otherwise unchanged
    """
    n = len(changestats_func_list)
    assert(len(theta) == n)
    accepted = 0
    addChangeStats = np.zeros(n)
    delChangeStats = np.zeros(n)
    for k in range(sampler_m):
        # basic sampler: select two nodes i and j uniformly at random
        # and toggle edge between them
        (i, j) = random.sample(list(range(G.numNodes())), 2) # without replacement
        assert(i != j)
        isDelete = (G.isArc(i, j))
        if isDelete:
            G.removeArc(i, j)

        # compute change statistics for each of the n statistics using the
        # list of change statistic functions
        changestats = np.zeros(n)
        for l in range(n):
            changestats[l] = changestats_func_list[l](G, i, j)
            assert(changestats[l] >= 0)
        changeSignMul = -1 if isDelete else +1
        total = np.sum(theta * changeSignMul * changestats)
        if random.uniform(0, 1) < math.exp(total):
            accepted += 1
            if performMove:
                # actually accept the move.
                # if deleting, we have already done it. For add move,
                # add the edge now
                if not isDelete:
                    G.insertArc(i, j)
            else:
                # if we are not to actually perform the moves, then reverse
                # changes for delete move made so G and twoPath same as before
                if isDelete:
                    G.insertArc(i, j)
            if isDelete:
                delChangeStats += changestats
            else:
                addChangeStats += changestats
        elif isDelete: # move not accepted, so reverse change for delete
            G.insertArc(i, j)

    acceptance_rate = float(accepted) / sampler_m
    return (acceptance_rate, addChangeStats, delChangeStats)


def algorithm_S(G, changestats_func_list, M1, theta_outfile):
    """

     Algorithm S

     Parameters:
        G                   - Digraph object for graph to estimate
        changestat_func_v   - vector of change statistics funcions
        M1                  - number of iterations of Algorithm S
        theta_outfile       - open for write file to write theta values

     Returns:
       tuple with:
         theta               - numpy vector of theta values at end
         Dmean               - derivative estimate value at end

    """
    ACA = 0.1 # multiplier of da to get K1A step size multiplier
    n = len(changestats_func_list)
    theta = np.zeros(n)
    D0 = np.zeros(n)
    for t in range(M1):
        accepted = 0
        (acceptance_rate,
         addChangeStats,
         delChangeStats) = BasicSampler(G, changestats_func_list, theta,
                                        performMove = False)
        dzA = delChangeStats - addChangeStats
        dzAmean = dzA / sampler_m
        sumChangeStats = addChangeStats + delChangeStats
        assert(np.all(sumChangeStats >= 0)) # zero is handled below
        D0 += dzA**2 # 1/D0 is squared derivative
        da = np.zeros(n)
        for l in range(n):
            if (sumChangeStats[l] != 0):
                da[l] = ACA  / sumChangeStats[l]**2
        theta_step = np.sign(dzAmean) * da * dzA**2
        MAXSTEP = 0.1 # limit maximum step size
        theta_step = np.where(theta_step > MAXSTEP, MAXSTEP, theta_step)
        theta += theta_step
        theta_outfile.write(str(t-M1) + ' ' + ' '.join([str(x) for x in theta])
                            + ' ' + str(acceptance_rate) + '\n')
    Dmean = sampler_m / D0
    return(theta, Dmean)
        

def algorithm_EE(G, changestats_func_list, theta, D0,
                 Mouter, M, theta_outfile, dzA_outfile):
    """
    Algorithm EE

    Parameters:
       G                   - Digraph object for graph to estimate
       changestats_func_list-list of change statistics funcions
       theta               - corresponding vector of initial theta values
       D0                  - corresponding vector of inital D0 from Algorithm S
       Mouter              - iterations of Algorithm EE (outer loop)
       M                   - iterations of Algorithm EE (inner loop)
       theta_outfile       - open for write file to write theta values
       dzA_outfile         - open for write file to write dzA values

     Returns:
         numpy vector of theta values at end

    """
    ACA = 1e-09    # multiplier of D0 to get step size multiplier da (K_A)
    compC = 1e-02  # multiplier of sd(theta)/mean(theta) to limit theta variance
    n = len(changestats_func_list)
    dzA = np.zeros(n)  # zero outside loop, dzA accumulates in loop
    t = 0
    for touter in range(Mouter):
        thetamatrix = np.empty((M, n)) # rows theta vectors, 1 per inner iter
        for tinner in range(M):
            accepted = 0
            (acceptance_rate,
             addChangeStats,
             delChangeStats) = BasicSampler(G, changestats_func_list, theta,
                                            performMove = True)
            dzA += addChangeStats - delChangeStats  # dzA accumulates here
            da = D0 * ACA
            theta_step = -np.sign(dzA) * da * dzA**2
            theta += theta_step
            theta_outfile.write(str(t) + ' ' + ' '.join([str(x) for x in theta]) + 
                                ' ' + str(acceptance_rate) + '\n')
            dzA_outfile.write(str(t) + ' ' + ' '.join([str(x) for x in dzA]) + '\n')
            thetamatrix[tinner,] = theta
            t += 1
        thetamean = np.mean(thetamatrix, axis = 0) # mean theta over inner loop
        thetasd   = np.std(thetamatrix, axis = 0)  # standard deviation
        thetamean = np.where(np.abs(thetamean) < 1, np.ones(n), thetamean) # enforce minimum magnitude 1 to stop sticking at zero
        DD = thetasd / np.abs(thetamean)
        D0 *= compC / DD # to limit size of fluctuations in theta (see S.I.)
        
    return theta

#-------------------------------------------------------------------------------


def run_on_network_attr(edgelist_filename, param_func_list, labels,
                        binattr_filename=None,
                        catattr_filename=None):
    """
    Run on specified network with binary and/or categorical attributes.
    
    Parameters:
         edgelist_filename - filename of Pajek format edgelist 
         param_func_list   - list of change statistic functions corresponding
                             to parameters to estimate
         labels            - list of strings corresponding to param_func_list
                             to label output (header line)
         binattr_filename - filename of binary attributes (node per line)
                            Default None, in which case no binary attr.
         catattr_filename - filename of categorical attributes (node per line)
                            Default None, in which case no categorical attr.
    Write output to ifd_theta_values_<basename>.txt and
                    ifd_dzA_values_<basename>.txt
    where <basename> is the baesname of edgelist filename e..g
    if edgelist_filename is edges.txt then ifd_theta_values_edges.txt
    and ifd_dzA_values_edges.txt
    WARNING: these files are overwritten.
    """
    assert(len(param_func_list) == len(labels))
    basename = os.path.splitext(os.path.basename(edgelist_filename))[0]
    THETA_OUTFILENAME = THETA_PREFIX + basename + os.extsep + 'txt'
    DZA_OUTFILENAME = DZA_PREFIX + basename  + os.extsep + 'txt'

    G = Digraph(edgelist_filename, binattr_filename, catattr_filename)

    M1_steps = 500
    # steps of Alg 1    
    M1 = int(M1_steps * G.density()*(1 - G.density())*G.numNodes()**2 / sampler_m)

    Mouter = 500 # outer iterations of Algorithm EE
    Msteps = 100 # multiplier for number of inner steps of Algorithm EE
    # inner steps of EE    
    M = int(Msteps * G.density()*(1 - G.density())*G.numNodes()**2 / sampler_m)

    print('M1 = ', M1, ' Mouter = ', Mouter, ' M = ', M)
    
    theta_outfile = open(THETA_OUTFILENAME, 'w',1) # 1 means line buffering
    theta_outfile.write('t ' + ' '.join(labels) + ' ' + 'AcceptanceRate' + '\n')
    print('Running Algorithm S...', end=' ')
    start = time.time()
    (theta, Dmean) = algorithm_S(G, param_func_list, M1, theta_outfile)
    print(time.time() - start, 's')
    print('after Algorithm S:')
    print('theta = ', theta)
    print('Dmean = ', Dmean)
    dzA_outfile = open(DZA_OUTFILENAME, 'w',1)
    dzA_outfile.write('t ' + ' '.join(labels) + '\n')
    print('Running Algorithm EE...', end=' ')
    start = time.time()
    theta = algorithm_EE(G, param_func_list, theta, Dmean,
                         Mouter, M, theta_outfile, dzA_outfile)
    print(time.time() - start, 's')
    theta_outfile.close()
    dzA_outfile.close()
    print('at end theta = ', theta)

    

def run_example():
    """
    example run on simulated 500 node network
    """
    run_on_network_attr(
        'sample_statistics_n500_directed_binattr_sim420000000.txt',
        [changeArc, changeReciprocity, changeAltInStars, changeAltOutStars,
         changeAltKTrianglesT, changeAltTwoPathsTD,
         changeReceiver, changeSender, changeInteraction],
        ["Arc", "Reciprocity", "AinS", "AoutS", "AT-T", "A2P-TD",
         "Receiver", "Sender", "Interaction"],
        'binaryAttribute_50_50_n500.txt')
