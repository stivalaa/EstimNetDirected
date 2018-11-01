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
## Usage: Rscript plotEsimtNetDirectedSimFit.R netfilename simNetFilePrefix
##  netfilename is the Pajek format observed graph (the input arclistFile
##     for EstimNetDirected)
##  simNetFilePreifx is the prefix of the simulated network filenames
##    this files have _x.net appended by EstimNetDirected, where x
##    is taks number.
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

library(igraph)

library(grid)
library(gridExtra)
library(ggplot2)
library(reshape)
library(doBy)
library(scales)


args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript plotEstimNetDirectedSimFit.R netfilename simNetFilePrefix\n")
  quit(save="no")
}
netfilename <- args[1]
simnetfileprefix <- args[2]

obscolour <- 'red' # colour to plot observed graph points/lines
## simulated graph statistics will be boxplot on same plot in default colour

graph_glob <- paste(simnetfileprefix, "_[0-9]*[.]net", sep='')
outfilename <- paste(simnetfileprefix, "pdf", sep='.')

g_obs <- read.graph(netfilename, format='pajek')

cat('Reading ', length(graph_glob), ' graphs...\n')
system.time(sim_graphs <- sapply(Sys.glob(graph_glob),
                                 FUN = function(f) read.graph(f, format='pajek'),
                                 simplify = FALSE))

components <- sapply(sim_graphs, function(g) length(decompose.graph(g)))
ccs <- sapply(sim_graphs, function(g) transitivity(g, type="global"))
    

plotlist <- list()

cat('obs components: ', length(decompose.graph(g_obs)), '\n')
cat('sim components: ', components, '\n')
p <- ggplot() + geom_boxplot(aes(x = 'components', y = components))
p <- p + geom_point(aes(x = as.numeric(ordered('components')),
                        y = length(decompose.graph(g_obs)),
                        colour = obscolour))
p <- p + theme(legend.position = 'none')
plotlist <- c(plotlist, list(p))

cc_obs <- transitivity(g_obs, type='global')
cat('obs transitivity: ', cc_obs, '\n')
cat('sim transitivity: ', ccs, '\n')
p <- ggplot() + geom_boxplot(aes(x = 'transitivity', y = ccs))
p <- p + geom_point(aes(x = as.numeric(ordered('transitivity')),
                        y = cc_obs,
                        colour = obscolour))
p <- p + ylab('global clustering coefficient')
p <- p + theme(legend.position = 'none')
plotlist <- c(plotlist, list(p))


pdf(outfilename, onefile=FALSE, paper="special", width=9, height=6)
do.call(grid.arrange, plotlist)
dev.off()
