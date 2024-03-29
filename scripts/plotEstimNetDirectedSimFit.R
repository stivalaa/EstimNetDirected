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
## Usage: Rscript plotEsimtNetDirectedSimFit.R [-s] [-c] [-g] [-d] [-t] netfilename simNetFilePrefix
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
##       and
##       <simNetFilePrefix>_cycledist.eps (if -y is also specified)
##
##  -c : use induced subgraphs as generated by extractcERGM2subgraphs.py
##       to do GoF test for citation ERGM (cERGM) as an altnerative to
##       plotEstimNetDirectedcERGMSimFit2.R. The reason to do this is
##       that R/igraph is far too slow to be practical for this task on
##       larger networks (e.g. taking several hours per simulated graph
##       to get the induced subgraphs for the MAREC patent data, versus
##       about 1 minute with Pyton and SNAP in extractcERGM2subgraphs.py
##
##  -g : do NOT do geodesic distance distribution plot (useful as this can
##       be unusuably slow, especially on large networks)
##
##  -d : do NOT do dyadwise shared partners plot (useful as even on bipartite
##       graph where edgewise shared partners not done [always zero], the dsp
##       can take vast amounts of memory (has to be killed on node with "only"
##       64 GB) so unusuably slow / too much resoureces.
##
##  -t : do bipartite clustering coefficients (computed with tnet
##       pacakge). Note this is not done by default as it is very slow.
##
##  -y MAX_CYCLELEN : do cycle length distrubition (up to MAX_CYCLELEN). Note this is
##       not done by default (value 0) as it can be very slow
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

library(optparse)

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

## for library(parallel) using simplify2array(mclapply(...)) instead of
## sapply(...) and mclapply(...) instead of lapply(...)
cat('mc.cores =', getOption("mc.cores"), '\n')

args <- commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("-s", "--subplots"), action="store_true", default=FALSE,
                 help="Also put subplots in separate files"),
  make_option(c("-c", "--cergm_mode"), action="store_true", default=FALSE,
                 help="cERGM mode, use induced subgraphs from extractcERGM2subgraphs.py"),
  make_option(c("-g", "--no_geodesics"), action="store_true", default=FALSE,
                 help="do not do geodesic distance distribution"),
  make_option(c("-d", "--no_dsp"), action="store_true", default=FALSE,
                 help="do not do dyadwise shared partner distribution"),
  make_option(c("-t", "--bipartite_cc"), action="store_true", default=FALSE,
                 help="do bipartite clustering coefficients"),
  make_option(c("-y", "--max_cyclelen"), type="integer", default=0,
                help="max cycle length [default %default]")
  )
parser <- OptionParser(usage = "%prog [options] netfilename simNetFilePrefix",
                       option_list = option_list)
arguments <- parse_args(parser, positional_arguments = 2)
opt <- arguments$options
args <- arguments$args

do_subplots <- opt$subplots
cergm_mode <- opt$cergm_mode
do_geodesic <- !opt$no_geodesics
do_dsp <- !opt$no_dsp
do_bipartite_cc <- opt$bipartite_cc
MAX_CYCLELEN <- opt$max_cyclelen
do_cycledist <- (MAX_CYCLELEN > 0)


netfilename <- args[1]
simnetfileprefix <- args[2]

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

if (is.bipartite(g_obs)) {
  cat("bipartite graph\n")
  stopifnot(all(sapply(sim_graphs, FUN = function(g) is.bipartite(g))))
}

if (cergm_mode) {
  cat("cERGM mode\n")
  #cat("observed number of nodes: ", vcount(g_obs), "\n")
  #cat("simulated number of nodes:", sapply(sim_graphs, function(g) vcount(g)), "\n")
} else {
  num_nodes <- vcount(g_obs)
  ## all simulated graphs must have the same number of nodes
  stopifnot(length(unique((sapply(sim_graphs, function(g) vcount(g))))) == 1)
  ## and it must be the same a the number of nodes in the observed graph
  stopifnot(num_nodes == vcount(sim_graphs[[1]]))
}

num_sim <- length(sim_graphs)


## build the list of plots
plotlist <- build_sim_fit_plots(g_obs, sim_graphs, do_subplots, do_geodesic,
                                do_dsp, do_bipartite_cc, do_cycledist,
                                MAX_CYCLELEN)

if (cergm_mode) {
  ## add plot for number of nodes in each induced subgraph, with
  ## boxplot for simulated graphs and red point for observed
  ## as receivers of an arc from the last time period, for each graph
  ##print(vcount(g_obs)) #XXX
  sim_nodecounts <- sapply(sim_graphs, vcount)
  ##summary(sim_nodecounts) #XXX
  xlabel <-'Nodes in induced subgraph'
  p <- ggplot() + geom_boxplot(aes(x = xlabel,
                                   y = sim_nodecounts))
  p <- p + geom_point(aes(x = as.numeric(ordered(xlabel)),
                          y = vcount(g_obs),
                          colour = obscolour))
  p <- p + ylab('Number of nodes')
  p <- p + ptheme + theme(axis.title.x = element_blank())
  plotlist <- c(plotlist, list(p))
}

###
### Write the plots to PDF
###
cat("writing plots to PDF file ", outfilename, "\n")
pdf(outfilename, onefile=FALSE, paper="special", width=9.9, height=6.6)
system.time(do.call(grid.arrange, plotlist))
dev.off()
