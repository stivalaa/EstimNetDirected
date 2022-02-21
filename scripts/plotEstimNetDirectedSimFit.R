#!/usr/bin/Rscript
##
## File:    plotEstimNetDirectedSimFit.R
## Author:  Alex Stivala
## Created: November 2018
##
## Plot some simple graph statistics of the observed network on same
## plot as those of the networks output from the end of the MCMC
## simulation (when the outputSimulatedNetworks = True parameter is
## set in the config file to write these to the files with prefix set
## by simNetFilePrefix in the config file).  This gives some
## preliminary "goodness-of-fit" check that the networks from the MCMC
## simulation process are consistent with the observed network, in
## addition to the dzA value plot from plotEstimNetDirectedResults.R
## The latter should be centred around zero showing the values of the
## statistics included in the model are not significantly different in
## the simulated networks, while this plot uses some statistics not in the
## model explicitly (mainly degree distribution) which for a good model
## should also not be significantly different from those of the observed
## network being modeled.
##
## This same script is used for standard goodness-of-fit plots for
## networks simulated from a model with the SimulateERGM program.
##
## Usage: Rscript plotEsimtNetDirectedSimFit.R [-s] netfilename simNetFilePrefix
##  netfilename is the Pajek format observed graph (the input arclistFile
##     for EstimNetDirected)
##  simNetFilePreifx is the prefix of the simulated network filenames
##    this files have _x.net appended by EstimNetDirected or SimulateERGM,
##    where is task or iteration number.
##
##  -s : do subplots separately. As well as pdf file, .eps file for certain
##       of the subplots e.g. triad census (log scale) is also done as 
##       separte files, with name of subplot appeneded e.g. as
##       <simNetFilePrefix>_triadcensus.eps
##
## Output file is simfitPrefix.pdf (where Prefix is the simNetFilePrefix).
## WARNING: output file is overwritten
##
## Example:
## Rscript plotEsimtNetDirectedSimFit.R ../pythonDemo/polblogs/polblogs_arclist.txt sim_polblogs
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

## Note edge-wise and dyad-wise shared partner distribution like
## statnet GoF or something similar cannot be done with R/igraph
## (although similarity.dice in Python/igraph could be useful as it
## has the pairs not just vids parameters, R/igraph does not). See:
## https://github.com/igraph/igraph/issues/331
## So therefore using statnet library to calculate this, so have
## to load intergraph
## http://mbojan.github.io/intergraph/
## to convert to Network object and statnet,
## also too slow to be used for larger networks
## Also note only loading statnet and intergraph if required (network is
## small enough that they are practical) is if they are loaded here
## it seems another R problem means we run out of memory even on a 64 GB
## limit even though not actually used and without them it worked in less
## than 8 GB. (Really should just do everything in Python again, wasting
## far too much time & effort with R being too slow and too many problems,
## ending up having to rewrite in Python anyway like for snowball sampling...)
##

library(igraph)

library(grid)
library(gridExtra)
library(ggplot2)
library(reshape)
library(doBy)
library(scales)


## read in R source file from directory where this script is located
##http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep=.Platform$file.sep))
}

source_local('snowballSample.R')
source_local('simFitPlots.R')




###
### Main
###

args <- commandArgs(trailingOnly=TRUE)
basearg <- 0
do_subplots <- FALSE
if (length(args) > 3) {
  cat("Usage: Rscript plotEstimNetDirectedSimFit.R [-s] netfilename simNetFilePrefix\n")
  quit(save="no")
} else if (length(args) == 3) {
  if (args[1] == "-s") {
    do_subplots <- TRUE
  } else {
  cat("Usage: Rscript plotEstimNetDirectedSimFit.R [-s] netfilename simNetFilePrefix\n")
  quit(save="no")
  }
  basearg <- basearg + 1
}
netfilename <- args[basearg+1]
simnetfileprefix <- args[basearg+2]


graph_glob <- paste(simnetfileprefix, "_[0-9]*[.]net", sep='')
outfilename <- paste(simnetfileprefix, "pdf", sep='.')

g_obs <- read_graph_file(netfilename)

sim_files <- Sys.glob(graph_glob)
if (length(sim_files) == 0) {
  ## try .gz (compressed) instead
  cat('No .net files with prefix ', simnetfileprefix, ' found, trying .net.gz instead\n')
  graph_glob <- paste(simnetfileprefix, "_[0-9]*[.]net.gz", sep='')
  sim_files <- Sys.glob(graph_glob)
}

cat('Reading ', length(sim_files), ' graphs...\n')
system.time(sim_graphs <- sapply(sim_files,
                                 FUN = function(f) read_graph_file(f),
                                 simplify = FALSE))

num_nodes <- vcount(g_obs)
## all simulated graphs must have the same number of nodes
stopifnot(length(unique((sapply(sim_graphs, function(g) vcount(g))))) == 1)
## and it must be the same a the number of nodes in the observed graph
stopifnot(num_nodes == vcount(sim_graphs[[1]]))
num_sim <- length(sim_graphs)

ptheme <-  theme(legend.position = 'none')

## build the list of plots
plotlist <- build_sim_fit_plots(g_obs, sim_graphs)

###
### Write the plots to PDF
###
cat("writing plots to PDF file ", outfilename, "\n")
pdf(outfilename, onefile=FALSE, paper="special", width=9.9, height=6.6)
system.time(do.call(grid.arrange, plotlist))
dev.off()
