#!/usr/bin/Rscript
##
## File:    statsEstimNetDirectedSimFit.R
## Author:  Alex Stivala
## Created: September 2021
##
## Compute some statistics related to the distance of the
## observed graphs statistics from those of simulated graphs.
##
## Usage: Rscript statsEsimtNetDirectedSimFit.R  netfilename simNetFilePrefix
##  netfilename is the Pajek format observed graph (the input arclistFile
##     for EstimNetDirected)
##  simNetFilePreifx is the prefix of the simulated network filenames
##    this files have _x.net appended by EstimNetDirected, where 
##    is task number.
##
## Example:
## Rscript statsEsimtNetDirectedSimFit.R ../pythonDemo/polblogs/polblogs_arclist.txt sim_polblogs
##  which will use input files sim_polblogs_0.net etc.
##
## Uses the igraph library to read Pajek format graph files and
## compute graph statistics:
##
##   Csardi G, Nepusz T: The igraph software package for complex network
##   research, InterJournal, Complex Systems
##   1695. 2006. http://igraph.org
##
##
## Note that EstimNetDirected can handle very large graphs, but R and igraph
## can be very slow and may not be able to practically work for larger
## graphs especially for triad census etc., but at least degree distribution
## should be able to be computed.
## (Also for speed and convenience this script reads all the graphs in
## so could use lots of memory).
##
##

library(igraph)
library(MASS) #for ginv()


## read in R source file from directory where this script is located
##http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep=.Platform$file.sep))
}

source_local('snowballSample.R')


##
## Return Mahalanobis distance of observed to simulated
## degree distribution, for in or out degree
##
## Parameters:
##    g_obs:       observed graph igraph object
##    sim_graphs:  simulated graphs list of igraph objects
##    mode:       'in' or 'out' for indegree or outdegree respectively
##
## Return value:
##    Mahalanobis distance between observed and simulated distributions
##
##
deg_distr_distance <- function(g_obs, sim_graphs, mode) {
    maxdeg <- max(sapply(sim_graphs, function(g) degree(g, mode=mode)),
                  degree(g_obs, mode=mode))
    deg_df <- data.frame(sim = rep(1:num_sim, each=(maxdeg+1)),
                           degree = rep(0:maxdeg, num_sim),
                           count = NA)
    # construct tmatrix where columns are counts for degree0, degree1, etc.
    # and rows are the observations (each row a simulated network)
    deg_matrix <- matrix(nrow = num_sim, ncol = (1+maxdeg))
    for (i in 1:num_sim) {
        ## https://stackoverflow.com/questions/1617061/include-levels-of-zero-count-in-result-of-table
        deg_table <- table(factor(degree(sim_graphs[[i]], mode = mode),
                                  levels=0:maxdeg))
        deg_matrix[i,] <- t(matrix(deg_table))
    }

    obs_deg_df <- data.frame(degree = rep(0:maxdeg),
                               count = NA)
    obs_deg_table <- table(factor(degree(g_obs, mode=mode), levels=0:maxdeg))
    ## convert to matrix where columns are degree0, degree1,... etc.
    ## and single row for the single observation
    obs_deg_matrix <- t(matrix(obs_deg_table))

    # Add observed data as last row
    deg_matrix <- rbind(deg_matrix, obs_deg_matrix)

    #print(deg_matrix)#XXX

    dcov <- cov(deg_matrix)
    ## http://r.789695.n4.nabble.com/Catching-errors-from-solve-with-near-singular-matrices-td4652794.html
    if (rcond(dcov) < .Machine$double.eps)  {
        cat("WARNING: ", mode,"-degree matrix is singular, using pseudo-inverse of covariance matrix\n", file=stderr())
        inverted_cov_deg_matrix <- ginv(dcov) # generalized (Moore-Penrose) inverse of dcov
    } else {
        inverted_cov_deg_matrix <- solve(dcov) # inverse of dcov
    }

    mdist <- mahalanobis(deg_matrix, colMeans(deg_matrix), inverted_cov_deg_matrix, inverted=TRUE)
    print(mdist)#XXX
    return(mdist[length(mdist)]) # Mahalanobis distance of observed (last in vector from mahalanobis())
}

##
## Return matrix where rows are
## degree distribution of each observation
## (simulated network, observed as last row), for in or out degree
##
## Parameters:
##    g_obs:       observed graph igraph object
##    sim_graphs:  simulated graphs list of igraph objects
##    maxdeg:     maximum degree
##    mode:       'in' or 'out' for indegree or outdegree respectively
##
## Return value:
##    Matrix of degree distributions
##
##
build_deg_matrix <- function(g_obs, sim_graphs, maxdeg, mode) {

    deg_df <- data.frame(sim = rep(1:num_sim, each=(maxdeg+1)),
                           degree = rep(0:maxdeg, num_sim),
                           count = NA)
    # construct matrix where columns are counts for degree0, degree1, etc.
    # and rows are the observations (each row a simulated network)
    deg_matrix <- matrix(nrow = num_sim, ncol = (1+maxdeg))
    for (i in 1:num_sim) {
        ## https://stackoverflow.com/questions/1617061/include-levels-of-zero-count-in-result-of-table
        deg_table <- table(factor(degree(sim_graphs[[i]], mode = mode),
                                  levels=0:maxdeg))
        deg_matrix[i,] <- t(matrix(deg_table))
    }

    obs_deg_df <- data.frame(degree = rep(0:maxdeg),
                               count = NA)
    obs_deg_table <- table(factor(degree(g_obs, mode=mode), levels=0:maxdeg))
    ## convert to matrix where columns are degree0, degree1,... etc.
    ## and single row for the single observation
    obs_deg_matrix <- t(matrix(obs_deg_table))

    # Add observed data as last row
    deg_matrix <- rbind(deg_matrix, obs_deg_matrix)

    return(deg_matrix)
}


###
### Main
###

args <- commandArgs(trailingOnly=TRUE)
basearg <- 0
if (length(args) != 2) {
  cat("Usage: Rscript statsEstimNetDirectedSimFit.R netfilename simNetFilePrefix\n")
  quit(save="no")
}
netfilename <- args[basearg+1]
simnetfileprefix <- args[basearg+2]


graph_glob <- paste(simnetfileprefix, "_[0-9]*[.]net", sep='')
outfilename <- paste(simnetfileprefix, "pdf", sep='.')

g_obs <- read_graph_file(netfilename, directed = TRUE)

sim_files <- Sys.glob(graph_glob)
sim_graphs <- sapply(sim_files,
                     FUN = function(f) read_graph_file(f,
                                                       directed=TRUE),
                     simplify = FALSE)

num_nodes <- vcount(g_obs)
## all simulated graphs must have the same number of nodes
stopifnot(length(unique((sapply(sim_graphs, function(g) vcount(g))))) == 1)
## and it must be the same a the number of nodes in the observed graph
stopifnot(num_nodes == vcount(sim_graphs[[1]]))
num_sim <- length(sim_graphs)

## do stats for degree 0, degree 1,..., maxdegree
maxdegree <- 5 # capping at this limit as if too many cols get singular matrix
#maxindegree  <- max(degree(g_obs, mode='in'))
#maxoutdegree <- max(degree(g_obs, mode='out'))
maxindegree <- maxdegree
maxoutdegree <- maxdegree

## In degree

stats_matrix <- build_deg_matrix(g_obs, sim_graphs, maxindegree, mode='in')

## Out degree

stats_matrix <- cbind(stats_matrix, build_deg_matrix(g_obs, sim_graphs, maxoutdegree, mode='out'))

print(stats_matrix)#XX

statscov <- cov(stats_matrix)
inverted_cov_stats_matrix <- solve(statscov) # inverse of statscov
mdist <- mahalanobis(stats_matrix, colMeans(stats_matrix), inverted_cov_stats_matrix, inverted=TRUE)

print(mdist)#XXX
obs_mdist <- mdist[length(mdist)] # Mahalanobis distance of observed (last in vector from mahalanobis())
cat("Mahalanobis = ", obs_mdist, "\n")
