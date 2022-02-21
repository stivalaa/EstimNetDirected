#!/usr/bin/Rscript
##
## File:    simFitPlots.R
## Author:  Alex Stivala
## Created: February 2022
##
## Functions to build goodness-of-fit plots for plotEstimNetDirectedSimFit.R
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

## some statistics are too slow to practically compute on large networks,
## these are just guesses (and certainly for geodesic a 1.6 million node
## network could not be computed in 4 hours for example).
MAX_SIZE_GEODESIC <- 100000 ## do not do shortest paths if more nodes than this
MAX_SIZE_ESP_DSP <-  100000 ## do not do shared partners if more nodes than this

# http://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
orig_scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}
my_scientific_10 <- function(x) {
# also remove + and leading 0 in exponennt
  parse( text=gsub("e", " %*% 10^", gsub("e[+]0", "e", scientific_format()(x))) )
   
}

##
## Return plot of degree distribution, for in or out degree
##
## Parameters:
##    g_obs:       observed graph igraph object
##    sim_graphs:  simulated graphs list of igraph objects
##    mode:       'in' or 'out' for indegree or outdegree respectively
##                or 'all' for undirected graph
##
## Return value:
##    ggplot2 object to add to plot list
##
##
deg_distr_plot <- function(g_obs, sim_graphs, mode) {
    start = Sys.time()
    maxdeg <- max(sapply(sim_graphs, function(g) degree(g, mode=mode)),
                  degree(g_obs, mode=mode))
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
    degreetype <- 'degree'
    if (mode == 'in' || mode =='out') {
      degreetype <- paste(mode, 'degree', sep='-')
    }
    p <- p + xlab(degreetype) + ylab('fraction of nodes')
    if (maxdeg > 200) {
        p <- p + scale_x_discrete(breaks = seq(0, maxdeg, by = 200))
    } else if (maxdeg > 50) {
        p <- p + scale_x_discrete(breaks = seq(0, maxdeg, by = 10))
    } else if (maxdeg > 20) {
        p <- p + scale_x_discrete(breaks = seq(0, maxdeg, by = 5))
    }
    end = Sys.time()
    cat(mode, "-degree plotting took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    return(p)
}

##
## Return histogram of degree distribution, for in or out degree
##
## Parameters:
##    g_obs:       observed graph igraph object
##    sim_graphs:  simulated graphs list of igraph objects
##    mode:       'in' or 'out' for indegree or outdegree respectively
##                 or all for undirected
##    use_log:    TRUE to do log degree
##
## Return value:
##    ggplot2 object to add to plot list
##
deg_hist_plot <- function(g_obs, sim_graphs, mode, use_log) {
    start <- Sys.time()
    if (use_log) {
        dobs <- data.frame(degree = log(degree(g_obs, mode=mode)),
                           group = 'obs')
    } else {
        dobs <- data.frame(degree = degree(g_obs, mode=mode),
                           group = 'obs')
    }
    ## get degrees of all simulated graphs in one histogram
    simdegrees <- as.vector(sapply(sim_graphs, function(g) degree(g, mode=mode)))
    if (use_log) {
        dsim <- data.frame(degree = log(simdegrees), group = 'sim')
    } else {
        dsim <- data.frame(degree = simdegrees, group = 'sim')
    }
    dat <- rbind(dobs, dsim)
    end <- Sys.time()
    cat(mode, "-degree histogram data frame construction took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    start <- Sys.time()
    ## https://stackoverflow.com/questions/29287614/r-normalize-then-plot-two-histograms-together-in-r
    p <- ggplot(dat, aes(degree, fill = group, colour = group)) +
        geom_histogram(aes(y = ..density..),
                       alpha = 0.4, position = 'identity', lwd = 0.2)
    degreetype <- 'degree'
    if (mode == 'in' || mode =='out') {
      degreetype <- paste(mode, 'degree', sep='-')
    }
    p <- p + xlab(paste(ifelse(use_log, "log ", ""), degreetype, sep=''))
    p <- p + theme(legend.title=element_blank(),
                   legend.position = c(0.9, 0.8))
    end <- Sys.time()
    cat(mode, "-degree histogram plotting took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    return(p)
}

