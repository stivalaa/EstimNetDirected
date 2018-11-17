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

## TODO edge-wise and dyad-wise shared partner distribution like
## statnet GoF or something similar. See:
## https://github.com/igraph/igraph/issues/331

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


##
## Return plot of degree distribution, for in or out degree
##
## Parameters:
##    g_obs:       observed graph igraph object
##    sim_graphs:  simulated graphs list of igraph objects
##    mode:       'in' or 'out' for indegree or outdegree respectively
##
## Return value:
##    ggplot2 object to add to plot list
##
## TODO work out how to make x axis labels better; if max degree is high
## it is just a mess, need to work out how to label every 10th or 100th
## tick mark or something.
## TODO also sometimes it is better to plot this with y on log scale as
## if there is something like a log-normal or power-law degree distribution
## the graph is not very useful to read without log transformation.
deg_distr_plot <- function(g_obs, sim_graphs, mode) {
    start = Sys.time()
    maxdeg <- max(sapply(sim_graphs, function(g) degree(g, mode=mode)),
                  degree(g_obs, mode='in'))
    cat("Max ", mode, " degree is ", maxdeg, "\n")
    deg_df <- data.frame(sim = rep(1:num_sim, each=(maxdeg+1)),
                           degree = rep(0:maxdeg, num_sim),
                           count = NA)
    end = Sys.time()
    cat(mode, "-degree init took ", as.numeric(difftime(end, start, unit="secs")),"s\n")
    start = Sys.time()
    for (i in 1:num_sim) {
        ## https://stackoverflow.com/questions/1617061/include-levels-of-zero-count-in-result-of-table
        deg_table <- table(factor(degree(sim_graphs[[i]], mode = mode),
                                  levels=0:maxdeg))
        deg_df[which(deg_df[,"sim"] == i), "count"] <- deg_table
    }
    deg_df$degree <- as.factor(deg_df$degree)
    deg_df$count[which(is.na(deg_df$count))] <- 0
    deg_df$nodefraction <- deg_df$count / num_nodes
    end = Sys.time()
    cat(mode, "-degree sim data frame construction took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    start = Sys.time()
    obs_deg_df <- data.frame(degree = rep(0:maxdeg),
                               count = NA)
    obs_deg_table <- table(factor(degree(g_obs, mode=mode), levels=0:maxdeg))
    obs_deg_df$count <- as.numeric(obs_deg_table)
    ## without as.numeric() above get error "Error: geom_line requires
    ## the following missing aesthetics: y" when the plot is finally
    ## printed at the end. Who knows why... even though printing the
    ## data frame and the computations below are apparently not
    ## affected by this at all (does not happen with the boxplot for
    ## simulated degree distribution)
    obs_deg_df$degree <- as.factor(obs_deg_df$degree)
    obs_deg_df$count[which(is.na(obs_deg_df$count))] <- 0
    obs_deg_df$nodefraction <- obs_deg_df$count / num_nodes
    ##print(obs_deg_df)#XXX
    end = Sys.time()
    cat(mode, "-degree obs data frame construction took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    start = Sys.time()
    p <- ggplot(deg_df, aes(x = degree, y = nodefraction)) + geom_boxplot()
    p <- p + geom_line(data = obs_deg_df, aes(x = degree, y = nodefraction,
                                              colour = obscolour,
                                              group = 1))
    ## the "group=1" is ncessary in the above line otherwise get error
    ## "geom_path: Each group consists of only one observation. Do you
    ## need to adjust the group aesthetic?" and it does not work.
    ## https://stackoverflow.com/questions/27082601/ggplot2-line-chart-gives-geom-path-each-group-consist-of-only-one-observation
    p <- p + ptheme
    p <- p + xlab(paste(mode, '-degree', sep='')) + ylab('fraction of nodes')
    end = Sys.time()
    cat(mode, "-degree plotting took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    return(p)
}


###
### Main
###

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

g_obs <- read_graph_file(netfilename, directed = TRUE)

sim_files <- Sys.glob(graph_glob)
cat('Reading ', length(sim_files), ' graphs...\n')
system.time(sim_graphs <- sapply(sim_files,
                                 FUN = function(f) read_graph_file(f,
                                                                   directed=TRUE),
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

system.time(plotlist <- c(plotlist,
                          list(deg_distr_plot(g_obs, sim_graphs, 'in'))))


###
### Out degree
###

system.time(plotlist <- c(plotlist,
                          list(deg_distr_plot(g_obs, sim_graphs, 'out'))))


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
### giant component size
###

system.time(giant_component_sizes <- sapply(sim_graphs,
                                           function(g) vcount(giant.component(g))))
giant_component_sizes <- giant_component_sizes / num_nodes
obs_gcsize <- vcount(giant.component(g_obs))
cat('obs giant component size: ', obs_gcsize, '\n')
cat('sim giant component size: ', giant_component_sizes, '\n')
p <- ggplot() + geom_boxplot(aes(x = 'giant component', y = giant_component_sizes))
p <- p + geom_point(aes(x = as.numeric(ordered('giant component')),
                        y = obs_gcsize / num_nodes,
                        colour = obscolour))
p <- p + ylab('fraction of nodes')
p <- p + ptheme +   theme(axis.title.x = element_blank())
plotlist <- c(plotlist, list(p))



###
### (weakly) Connected components size distribution
###
## This is commented out as it always is too hard to read, a histogram
## works better but hard to superimpose these with multiple simluated graphs

## sim_component_sizes <- unlist(sapply(sim_graphs,
##                               function(g) sapply(decompose.graph(g),
##                                                  function(h) vcount(h))))
## obs_component_sizes <- sapply(decompose.graph(g_obs), function(g) vcount(g))
## maxcomponentsize <- max(c(sim_component_sizes, obs_component_sizes))
## componentsize_df <- data.frame(sim = rep(1:num_sim, each=maxcomponentsize),
##                                componentsize = rep(1:maxcomponentsize, num_sim),
##                                count = NA)
## for (i in 1:num_sim) {
##     ## using ineffecient and inelegant double loops as per above for degree
##     ## distribution. (see also necessity for as.character(j) below)
##     sim_component_sizes <- sapply(decompose.graph(sim_graphs[[i]]),
##                                   function(g) vcount(g))
##     componentsize_table <- table(sim_component_sizes)
##     for (j in 0:maxcomponentsize) {
##         componentsize_df[which(componentsize_df[,"sim"] == i &
##                                componentsize_df[,"componentsize"] == j,
##                                arr.ind=TRUE),
##                          "count"] <-
##             componentsize_table[as.character(j)]
##     }
## }
## componentsize_df$componentsize <- as.factor(componentsize_df$componentsize)
## componentsize_df$count[which(is.na(componentsize_df$count))] <- 0
## componentsize_df$nodefraction <- componentsize_df$count / num_nodes
## obs_componentsize_df <- data.frame(componentsize = 1:maxcomponentsize,
##                                    count = NA)
## obs_componentsize_table <- table(obs_component_sizes)
## for (j in 0:maxcomponentsize) {
##     obs_componentsize_df[which(obs_componentsize_df[,"componentsize"] == j,
##                                arr.ind=TRUE), "count"] <-
##         componentsize_table[as.character(j)]
## }
## obs_componentsize_df$componentsize <- as.factor(obs_componentsize_df$componentsize)
## obs_componentsize_df$count[which(is.na(obs_componentsize_df$count))] <- 0
## obs_componentsize_df$nodefraction <- obs_componentsize_df$count / num_nodes
## p <- ggplot(componentsize_df, aes(componentsize, nodefraction)) + geom_boxplot()
## p <- p + geom_line(data = obs_componentsize_df, aes(componentsize, nodefraction,
##                                             colour = obscolour,
##                                             group = 1))
## p <- p + ptheme
## p <- p + xlab('component size') + ylab('fraction of nodes')
## ##p <- p + scale_y_log10()
## plotlist <- c(plotlist, list(p))



###
### Transitivity (global clustering coefficient and avg. local clustering coef.)
###

cctypes <- c('global', 'average local')
system.time(ccs <- sapply(sim_graphs, function(g) transitivity(g, type="global")))
system.time(cc_obs <- transitivity(g_obs, type='global'))
cat('obs global cc: ', cc_obs, '\n')
cat('sim global cc: ', ccs, '\n')
system.time(ccs_localavg <- sapply(sim_graphs, function(g)
    transitivity(g, type='localaverage')))
system.time(cc_localavg_obs <- transitivity(g_obs, type='localaverage'))
p <- ggplot() + geom_boxplot(aes(x = factor('global', levels=cctypes), y = ccs))
p <- p + geom_point(aes(x = as.numeric(factor('global', levels=cctypes)),
                        y = cc_obs,
                        colour = obscolour))
p <- p + geom_boxplot(aes(x = factor('average local', levels=cctypes),
                          y = ccs_localavg))
p <- p + geom_point(aes(x = as.numeric(factor('average local', levels=cctypes)),
                        y = cc_localavg_obs,
                        colour = obscolour))
p <- p + ylab('clustering coefficient') + ptheme +
  theme(axis.title.x = element_blank())
p <- p + ylim(0, 1)
plotlist <- c(plotlist, list(p))


###
### Triad census
###

## Note that on large networks, triad census counts can overflow and
## give negative numbers
## https://github.com/igraph/igraph/issues/625
## https://github.com/igraph/igraph/issues/497

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
cat("writing plots to PDF file ", outfilename, "\n")
pdf(outfilename, onefile=FALSE, paper="special", width=9, height=6)
system.time(do.call(grid.arrange, plotlist))
dev.off()
