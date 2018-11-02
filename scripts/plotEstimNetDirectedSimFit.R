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

ptheme <-  theme(legend.position = 'none')


plotlist <- list()

###
### In degree
###
maxindeg <- max(sapply(sim_graphs, function(g) degree(g, mode='in')))
indeg_df <- data.frame(sim = rep(1:num_sim, each=(maxindeg+1)),
                       indegree = rep(0:maxindeg, num_sim),
                       count = NA)
for (i in 1:num_sim) {
    ## using inefficient and inelegant double loops as could not get
    ## replacement of all indegree values (for sim == i) of data frame
    ## to work, always get error "replacement has x rows, data has y"
    ## where y is total rows in data frame, not the subset, even
    ## though printing nrow showed correct z < y rows. Just too much
    ## time wasted trying to work out errors in R, gave up and did it
    ## this way.
    indeg_table <- table(degree(sim_graphs[[i]], mode = 'in'))
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
obs_indeg_df <- data.frame(indegree = rep(0:maxindeg),
                           count = NA)
obs_indeg_table <- table(degree(g_obs, mode='in'))
for (j in 0:maxindeg) {
    obs_indeg_df[which(obs_indeg_df[,"indegree"] == j, arr.ind=TRUE), "count"] <-
        indeg_table[as.character(j)]
}
obs_indeg_df$indegree <- as.factor(obs_indeg_df$indegree)
obs_indeg_df$nodefraction <- obs_indeg_df$count / num_nodes
p <- ggplot(indeg_df, aes(indegree, nodefraction)) + geom_boxplot()
p <- p + geom_line(data = obs_indeg_df, aes(indegree, nodefraction,
                                            colour = obscolour,
                                            group = 1))
## the "group=1" is ncessary in the above line otherwise get error
## "geom_path: Each group consists of only one observation. Do you
## need to adjust the group aesthetic?" and it does not work.
## https://stackoverflow.com/questions/27082601/ggplot2-line-chart-gives-geom-path-each-group-consist-of-only-one-observation

p <- p + ptheme
p <- p + xlab('in-degree') + ylab('fraction of nodes')
plotlist <- c(plotlist, list(p))





###
### Out degree
###
maxoutdeg <- max(sapply(sim_graphs, function(g) degree(g, mode='out')))
outdeg_df <- data.frame(sim = rep(1:num_sim, each=(maxoutdeg+1)),
                       outdegree = rep(0:maxoutdeg, num_sim),
                       count = NA)
for (i in 1:num_sim) {
    ## using inefficient and inelegant double loops as could not get
    ## replacement of all outdegree values (for sim == i) of data frame
    ## to work, always get error "replacement has x rows, data has y"
    ## where y is total rows in data frame, not the subset, even
    ## though printing nrow showed correct z < y rows. Just too much
    ## time wasted trying to work out errors in R, gave up and did it
    ## this way.
    outdeg_table <- table(degree(sim_graphs[[i]], mode = 'out'))
    for (j in 0:maxoutdeg) {
        outdeg_df[which(outdeg_df[,"sim"] == i &
                       outdeg_df[,"outdegree"] == j, arr.ind=TRUE), "count"] <-
            outdeg_table[as.character(j)]
        ## NB absolutely necessary to use as.character(j) in the line above
        ## otherwise it appears to work and has no errors/warnings but is wrong
        ## https://www.r-bloggers.com/indexing-with-factors/
    }
}
outdeg_df$outdegree <- as.factor(outdeg_df$outdegree)
outdeg_df$nodefraction <- outdeg_df$count / num_nodes
obs_outdeg_df <- data.frame(outdegree = rep(0:maxoutdeg),
                           count = NA)
obs_outdeg_table <- table(degree(g_obs, mode='out'))
for (j in 0:maxoutdeg) {
    obs_outdeg_df[which(obs_outdeg_df[,"outdegree"] == j, arr.ind=TRUE), "count"] <-
        outdeg_table[as.character(j)]
}
obs_outdeg_df$outdegree <- as.factor(obs_outdeg_df$outdegree)
obs_outdeg_df$nodefraction <- obs_outdeg_df$count / num_nodes
p <- ggplot(outdeg_df, aes(outdegree, nodefraction)) + geom_boxplot()
p <- p + geom_line(data = obs_outdeg_df, aes(outdegree, nodefraction,
                                            colour = obscolour,
                                            group = 1))
## the "group=1" is ncessary in the above line otherwise get error
## "geom_path: Each group consists of only one observation. Do you
## need to adjust the group aesthetic?" and it does not work.
## https://stackoverflow.com/questions/27082601/ggplot2-line-chart-gives-geom-path-each-group-consist-of-only-one-observation

p <- p + ptheme
p <- p + xlab('out-degree') + ylab('fraction of nodes')
plotlist <- c(plotlist, list(p))


###
### (weakly) Connected components
###

system.time(components <- sapply(sim_graphs, function(g) length(decompose.graph(g))))


cat('obs components: ', length(decompose.graph(g_obs)), '\n')
cat('sim components: ', components, '\n')
p <- ggplot() + geom_boxplot(aes(x = 'components', y = components))
p <- p + geom_point(aes(x = as.numeric(ordered('components')),
                        y = length(decompose.graph(g_obs)),
                        colour = obscolour))
p <- p + ptheme +   theme(axis.title.x = element_blank())
plotlist <- c(plotlist, list(p))


###
### Transitivity (global clustering coefficient)
###

system.time(ccs <- sapply(sim_graphs, function(g) transitivity(g, type="global")))
system.time(cc_obs <- transitivity(g_obs, type='global'))
cat('obs transitivity: ', cc_obs, '\n')
cat('sim transitivity: ', ccs, '\n')
p <- ggplot() + geom_boxplot(aes(x = 'transitivity', y = ccs))
p <- p + geom_point(aes(x = as.numeric(ordered('transitivity')),
                        y = cc_obs,
                        colour = obscolour))
p <- p + ylab('global clustering coefficient') + ptheme +
    theme(axis.title.x = element_blank())
plotlist <- c(plotlist, list(p))


###
### Triad census
###

nTriads <- choose(num_nodes, 3)
system.time(obs_triadcensus <- triad.census(g_obs))
num_triad_types <- length(obs_triadcensus)
stopifnot(num_triad_types == 16)
triadnames <- c('003', '012', '102', '021D', '021U', '021C', '111D',
                '111U', '030T', '030C', '201', '120D', '120U', '120C',
                '210', '300')
stopifnot(length(triadnames) == num_triad_types)
names(obs_triadcensus) <- triadnames
cat('obs triad census: ', obs_triadcensus, '\n')
sim_triadcensus_df <- data.frame(sim = rep(1:num_sim, each = num_triad_types),
                                 triad = rep(triadnames, num_sim),
                                 count = NA)
obs_triadcensus_df <- data.frame(triad = triadnames,
                                 count = NA)
## as for degree distributions, using loops as trying to do it "properly"
## in R was just too difficult
for (tname in triadnames) {
    obs_triadcensus_df[which(obs_triadcensus_df[,"triad"] == tname,
                             arr.ind=TRUE), "count"] <-
        obs_triadcensus[tname]
}
obs_triadcensus_df$triadfraction <- obs_triadcensus_df$count / nTriads
for (i in 1:num_sim) {
    system.time(sim_triadcensus <- triad_census(sim_graphs[[i]]))
    names(sim_triadcensus) <- triadnames
    cat('sim triad census ', i, ': ', sim_triadcensus, '\n')
    for (tname in triadnames) {
        sim_triadcensus_df[which(sim_triadcensus_df[,"sim"] == i &
                                 sim_triadcensus_df[,"triad"] == tname,
                                 arr.ind=TRUE), "count"] <-
            sim_triadcensus[tname]
    }
}

## make factor with triad names explicitly specified to keep them in order
sim_triadcensus_df$triad <- factor(sim_triadcensus_df$triad, levels = triadnames)
obs_triadcensus_df$triad <- factor(obs_triadcensus_df$triad, levels = triadnames)

sim_triadcensus_df$triadfraction <- sim_triadcensus_df$count / nTriads
## Remove 003 triad (empty graph) as it has has far larger fractoin than
## others so makes graph too hard to read (like in statnet GoF plots,
## everythign else is squashed too close to zero in comparison)
## sim_triadcensus_df <- sim_triadcensus_df[which(sim_triadcensus_df$triad != "003"),]
## obs_triadcensus_df <- obs_triadcensus_df[which(obs_triadcensus_df$triad != "003"),]
p <- ggplot(sim_triadcensus_df, aes(x = triad, y = triadfraction))
p <- p + geom_boxplot()
p <- p + ylab('fraction of traids') + ptheme +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
## who knows why the hjust and vjust are needed, or what they should
## be, but they do seem to be, otherwise labels are not positioned right
## (note depends on which versoin of R/ggplot2 being used, but this worked
## when I wrote it with R 3.4.2 ggplot2 2.2.1 on Windows 10 cygwin:
## https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
p <- p + geom_line(data = obs_triadcensus_df, aes(x = triad, y = triadfraction,
                                                  colour = obscolour,
                                                  group = 1))
p <- p + scale_y_log10()
plotlist <- c(plotlist, list(p))



###
### Write the plot to PDF
###

pdf(outfilename, onefile=FALSE, paper="special", width=9, height=6)
do.call(grid.arrange, plotlist)
dev.off()
