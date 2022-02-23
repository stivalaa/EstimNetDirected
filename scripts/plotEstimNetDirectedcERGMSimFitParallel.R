#!/usr/bin/Rscript
##
## File:    plotEstimNetDirectedcERGMSimFitParallel.R
## Author:  Alex Stivala
## Created: February 2022
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
## This version is specifically for the "citation ERGM" (cERGM)
## model (Schmid, Chen, & Desmarais 2021). Although in that paper
## the conventional statnet GoF is used, this is a new alternative method.
## It involves taking the subgraph induced by the nodes in the last
## time period, plus the union of all nodes that receive arcs from from
## them over all simulated networks (plus the observed), and only
## looking at the graph statistics of those networks, so they 
## aren't just hidden in the statistics of the entire graph.
##
## Schmid, C., Chen, T., & Desmarais, B. (2021). Generative Dynamics of 
## Supreme Court Citations: Analysis with a New Statistical Network Model. 
## Political Analysis, 1-20. doi:10.1017/pan.2021.20
##
##
## This version copied from plotEstimNetDirectedcERGMSimFit.R
## but using mclapply() in parallel library replacingly lapply() 
## of igraph neighbors() function over nodes in last term, as it
## is extrememely slow (as in about 4 hours for approx. 100 000 nodes).
##
## Usage: Rscript plotEsimtNetDirectedcERGMSimFitParallel.R [-s] netfilename termfilename simNetFilePrefix
##  netfilename is the Pajek format observed graph (the input arclistFile
##     for EstimNetDirected)
##  termfilename is the time periods (terms) of the nodes in the same
##     format as any node attribute for EstimNetDirected. I.e. the first
##     line is the attribute name (must be 'term' for this one) and each
##     subsequent line is the value for a node, in order of the node
##     numbering used in the arc list file (netfilename) i.e. first line
##     after header has term (time period) for node 0, the next line for node
##     1, and so on. The terms themsevles are numbered from 0 for the first
##     term (time period). This is the same format used for terms in cERGM
##     estimation in EstimNetDirected (and simulation in SimulateERGM); see
##     add_Cergm_terms_to_digraph() in src/graph.c
##  simNetFilePreifx is the prefix of the simulated network filenames
##    this files have _x.net appended by EstimNetDirected or SimulateERGM,
##    where is task or iteration number.
##
##  -s : do subplots separately. As well as pdf file, .eps file for certain
##       of the subplots e.g. triad census (log scale) is also done as 
##       separte files, with name of subplot appeneded e.g. as
##       <simNetFilePrefix>_triadcensus.eps
##
## Output file is simfitPrefix_cergm2.pdf (where Prefix is the simNetFilePrefix).
## WARNING: output file is overwritten
##
## Example:
## Rscript plotEsimtNetDirectedcERGMSimFitParallel.R ../pythonDemo/polblogs/polblogs_arclist.txt sim_polblogs
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

library(parallel)


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

cat('mc.cores =', getOption("mc.cores"), '\n')

args <- commandArgs(trailingOnly=TRUE)
basearg <- 0
do_subplots <- FALSE
if (length(args) > 4 || length(args) < 3) {
  cat("Usage: Rscript plotEstimNetDirectedcERGMSimFitParallel.R [-s] netfilename termfilename simNetFilePrefix\n")
  quit(save="no")
} else if (length(args) == 4) {
  if (args[1] == "-s") {
    do_subplots <- TRUE
  } else {
  cat("Usage: Rscript plotEstimNetDirectedcERGMSimFitParallel.R [-s] netfilename termfilename simNetFilePrefix\n")
  quit(save="no")
  }
  basearg <- basearg + 1
}
netfilename <- args[basearg+1]
termfilename <- args[basearg+2]
simnetfileprefix <- args[basearg+3]



graph_glob <- paste(simnetfileprefix, "_[0-9]*[.]net", sep='')
outfilename <- paste(paste(simnetfileprefix, 'cergm2', sep='_'), "pdf", sep='.')

g_obs <- read_graph_file(netfilename)
term_df <- read.table(termfilename, header=TRUE, stringsAsFactors=FALSE)
stopifnot(vcount(g_obs) == nrow(term_df))
stopifnot(min(term_df$term) == 0)
stopifnot(length(unique(term_df$term)) == 1 + max(term_df$term))
V(g_obs)$name <- 1:vcount(g_obs) # make sure graph is named so not relying on internal node ids
V(g_obs)$term <- term_df$term

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

## add term as node attriute to all simulated graphs
## also make sure graphs are named so not relying on internal node ids
for (i in 1:length(sim_graphs)) {
  V(sim_graphs[[i]])$term <- term_df$term
  V(sim_graphs[[i]])$name <- 1:vcount(sim_graphs[[i]])
}

## Get subgraph induced by union of nodes in last term and 
## all nodes (in all simulated and observed) that
## receive arcs from nodes in last term
maxterm <- max(term_df$term)
cat('maxterm = ', maxterm, '\n')
maxterm_nodes <- V(g_obs)[which(V(g_obs)$term == maxterm)]
cat('There are ', length(maxterm_nodes), ' nodes in last time period\n')
## neighbors() only takes a single node not a node sequence, so have to use
## Reduce(union, lapply(nodesequence), ...)
print(system.time(
maxterm_receiver_nodes <- Reduce(union, mclapply(maxterm_nodes, function(v)
                                    Filter(function(x) !(x %in% maxterm_nodes),
                                              neighbors(g_obs, v, mode='out'))))
))
warnings()#XXX
cat('There are ', length(maxterm_receiver_nodes), ' additional nodes in g_obs receiving arcs from nodes in last time period\n')
obs_extra_nodes <- length(maxterm_receiver_nodes)
## convert vertex sequences to ordinary vectors of integers (node name attr)
## so that they are compatible between different graph objects
maxterm_nodes <- as_ids(maxterm_nodes)
maxterm_receiver_nodes <- as_ids(maxterm_receiver_nodes)
sim_extra_nodes <- rep(NA, length(sim_graphs))
for (i in 1:length(sim_graphs)) {
  this_maxterm_nodes <- V(sim_graphs[[i]])[which(V(sim_graphs[[i]])$term == maxterm)]
  stopifnot( length(maxterm_nodes) == length(this_maxterm_nodes) && all(as_ids(this_maxterm_nodes) == maxterm_nodes) )
print(system.time(
  this_maxterm_receiver_nodes <- Reduce(union, mclapply(this_maxterm_nodes,
                                function(v)
                                Filter(function(x) !(x %in% this_maxterm_nodes),
                                    neighbors(sim_graphs[[i]], v, mode='out'))))
))
  warnings()#XXX
  cat('There are ', length(this_maxterm_receiver_nodes), ' additional nodes in simulated graph ', i, ' receiving arcs from nodes in last time period\n')
  sim_extra_nodes[i] <- length(this_maxterm_receiver_nodes)
  maxterm_receiver_nodes <- base::union(maxterm_receiver_nodes,
                                        as_ids(this_maxterm_receiver_nodes))
  cat('There are now ', length(maxterm_receiver_nodes), ' total unique extra nodes receiving ties from last time period\n')
}

g_obs <- induced.subgraph(g_obs, union(maxterm_nodes, maxterm_receiver_nodes))
cat('Final induced subgraphs have ', vcount(g_obs), ' nodes\n')
for (i in 1:length(sim_graphs)) {
  sim_graphs[[i]] <- induced.subgraph(sim_graphs[[i]],
                                 union(maxterm_nodes, maxterm_receiver_nodes))
}

## build the list of plots
plotlist <- build_sim_fit_plots(g_obs, sim_graphs)

## add plot for number of extra nodes in earlier time periods added
## as receivers of an arc from the last time period, for each graph
summary(sim_extra_nodes)#XXX
p <- ggplot() + geom_boxplot(aes(x = 'cited earlier term nodes',
                                 y = sim_extra_nodes))
p <- p + geom_point(aes(x = as.numeric(ordered('cited earlier term nodes')),
                        y = obs_extra_nodes,
                        colour = obscolour))
p <- p + ylab('number of nodes')
p <- p + ptheme + theme(axis.title.x = element_blank())
plotlist <- c(plotlist, list(p))

###
### Write the plots to PDF
###
cat("writing plots to PDF file ", outfilename, "\n")
pdf(outfilename, onefile=FALSE, paper="special", width=9.9, height=6.6)
system.time(do.call(grid.arrange, plotlist))
dev.off()
