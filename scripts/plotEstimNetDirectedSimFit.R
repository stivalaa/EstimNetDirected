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

sim_files <- Sys.glob(graph_glob)
cat('Reading ', length(sim_files), ' graphs...\n')
system.time(sim_graphs <- sapply(sim_files,
                                 FUN = function(f) read.graph(f, format='pajek'),
                                 simplify = FALSE))

num_nodes <- vcount(g_obs)
## all simulated graphs must have the same number of nodes
stopifnot(length(unique((sapply(sim_graphs, function(g) vcount(g))))) == 1)
## and it must be the same a the number of nodes in the observed graph
stopifnot(num_nodes == vcount(sim_graphs[[1]]))
num_sim <- length(sim_graphs)

ptheme <-  theme(legend.position = 'none',
                 axis.title.x = element_blank())

plotlist <- list()

indegree_obs <- table(degree(g_obs, mode='in'))
indegree_obs <- indegree_obs / num_nodes
##indegree_sim <- table((sapply(sim_graphs, function(g) degree(g, mode='in'))))
##indegree_sim <- indegree_sim / (num_nodes * num_sim)
maxindeg <- max(sapply(sim_graphs, function(g) degree(g, mode='in')))
indeg_df <- data.frame(sim = rep(1:num_sim, each=(maxindeg+1)),
                       indegree = rep(0:maxindeg, num_sim),
                       count = rep(NA, num_sim))
for (i in 1:num_sim) {
    ## using inefficient and inelegant double loops as could not get
    ## replacement of all indegree values (for sim == i) of data frame
    ## to work, always get error "replacement has x rows, data has y"
    ## where y is total rows in data frame, not the subset, even
    ## though printing nrow showed correct z < y rows. Just too much
    ## time wasted trying to work out errors in R, gave up and did it
    ## this way.
    indeg_table <- table(degree(sim_graphs[[i]], mode = 'in'))
    print(indeg_table) #XXX
    for (j in 0:maxindeg) {
        indeg_df[which(indeg_df[,"sim"] == i &
                       indeg_df[,"indegree"] == j, arr.ind=TRUE), "count"] <-
            indeg_table[as.character(j)]
        ## NB absolutely necessary to use as.character(j) in the line above
        ## otherwise it appears to work and has no errors/warnings but is wrong
        ## https://www.r-bloggers.com/indexing-with-factors/
    }
}
indeg_df$indegree <- as.factor(indeg_df$indegree)
indeg_df$nodefraction <- indeg_df$count / num_nodes
print(indeg_df)#XXX
p <- ggplot(indeg_df, aes(indegree, nodefraction)) + geom_boxplot()
## p <- p + geom_line(aes(x = as.numeric(ordered(names(indegree_obs))),
##                        y = as.numeric(indegree_obs), colour = obscolour))
p <- p + ptheme
p <- p + xlab('in-degree') + ylab('fraction of nodes')
plotlist <- c(plotlist, list(p))


system.time(components <- sapply(sim_graphs, function(g) length(decompose.graph(g))))
system.time(ccs <- sapply(sim_graphs, function(g) transitivity(g, type="global")))


cat('obs components: ', length(decompose.graph(g_obs)), '\n')
cat('sim components: ', components, '\n')
p <- ggplot() + geom_boxplot(aes(x = 'components', y = components))
p <- p + geom_point(aes(x = as.numeric(ordered('components')),
                        y = length(decompose.graph(g_obs)),
                        colour = obscolour))
p <- p + ptheme
plotlist <- c(plotlist, list(p))

cc_obs <- transitivity(g_obs, type='global')
cat('obs transitivity: ', cc_obs, '\n')
cat('sim transitivity: ', ccs, '\n')
p <- ggplot() + geom_boxplot(aes(x = 'transitivity', y = ccs))
p <- p + geom_point(aes(x = as.numeric(ordered('transitivity')),
                        y = cc_obs,
                        colour = obscolour))
p <- p + ylab('global clustering coefficient') + ptheme
plotlist <- c(plotlist, list(p))


pdf(outfilename, onefile=FALSE, paper="special", width=9, height=6)
do.call(grid.arrange, plotlist)
dev.off()
