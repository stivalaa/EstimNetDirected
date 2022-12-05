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

## for using simplify2array(mclapply(...)) instead of
## sapply(...) and mclapply(...) instead of lapply(...)
library(parallel)

## some statistics are too slow to practically compute on large networks,
## these are just guesses (and certainly for geodesic a 1.6 million node
## network could not be computed in 4 hours for example).
## TODO these are arbitrary, and the slowness has as much to do with
## density and network structure, not just number of nodes, so should find
## some better heuristic for this.
MAX_SIZE_GEODESIC <- 500000 ## do not do shortest paths if more nodes than this
MAX_SIZE_ESP_DSP <- 1000000 ## do not do shared partners if more nodes than this

## for cycle length distribution (only used if do_cycledist = True)
MAX_CYCLELEN <- 8

obscolour <- 'red' # colour to plot observed graph points/lines
## simulated graph statistics will be boxplot on same plot in default colour

ptheme <-  theme(legend.position = 'none')

# http://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
orig_scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}
my_scientific_10 <- function(x) {
# also remove + and leading 0 in exponennt
  parse( text=gsub("e", " %*% 10^", gsub("e[+]0", "e", scientific_format()(x))) )
   
}


##
## giant_component() - return largest connected component of the graph
## 
## Paramters:
##    graph - igraph to get giant componetn of
##
## Return value:
##    largest connected component of graph
##
giant.component <- function(graph) {
  cl <- clusters(graph)
  return(induced.subgraph(graph, which(cl$membership == which.max(cl$csize))))
}



##
## Return plot of degree distribution, for in or out degree
##
## Parameters:
##    g_obs:       observed graph igraph object
##    sim_graphs:  simulated graphs list of igraph objects
##    mode:       'in' or 'out' for indegree or outdegree respectively
##                or 'all' for undirected graph
##    btype:      igraph bipartite node type FALSE or TRUE, or NULL
##                if not bipartite. (Default NULL)
##
## Return value:
##    ggplot2 object to add to plot list
##
##
deg_distr_plot <- function(g_obs, sim_graphs, mode, btype=NULL) {
    num_sim <- length(sim_graphs)
    start = Sys.time()
    if (is.bipartite(g_obs)) {
      maxdeg <- max(unlist(sapply(sim_graphs,
           function(g) degree(g, V(g)[which(V(g)$type == btype)], mode=mode))),
           degree(g_obs, V(g_obs)[which(V(g_obs)$type == btype)], mode=mode))
    } else {
      maxdeg <- max(unlist(sapply(sim_graphs,
                                  function(g) degree(g, mode=mode))),
                    degree(g_obs, mode=mode))
    }
    cat("Max ", mode, " degree is ", maxdeg, "\n")
    deg_df <- data.frame(sim = rep(1:num_sim, each=(maxdeg+1)),
                           degree = rep(0:maxdeg, num_sim),
                           count = NA)
    end = Sys.time()
    cat(mode, "-degree init took ", as.numeric(difftime(end, start, unit="secs")),"s\n")
    start = Sys.time()
    for (i in 1:num_sim) {
      ## https://stackoverflow.com/questions/1617061/include-levels-of-zero-count-in-result-of-table
      if (is.bipartite(g_obs)) {
        deg_table <- table(factor(degree(sim_graphs[[i]],
    V(sim_graphs[[i]])[which(V(sim_graphs[[i]])$type == btype)], mode = mode),
                                  levels=0:maxdeg))
      } else {
        deg_table <- table(factor(degree(sim_graphs[[i]], mode = mode),
                                  levels=0:maxdeg))
      }
      deg_df[which(deg_df[,"sim"] == i), "count"] <- deg_table
    }
    deg_df$degree <- as.factor(deg_df$degree)
    deg_df$count[which(is.na(deg_df$count))] <- 0
    deg_df$nodefraction <- deg_df$count / sapply(sim_graphs, vcount)
    end = Sys.time()
    cat(mode, "-degree sim data frame construction took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    start = Sys.time()
    obs_deg_df <- data.frame(degree = rep(0:maxdeg),
                               count = NA)
    if (is.bipartite(g_obs)) {
      obs_deg_table <- table(factor(degree(g_obs,
                          V(g_obs)[which(V(g_obs)$type == btype)], mode=mode),
         levels=0:maxdeg))
    } else {
      obs_deg_table <- table(factor(degree(g_obs, mode=mode), levels=0:maxdeg))
    }
    obs_deg_df$count <- as.numeric(obs_deg_table)
    ## without as.numeric() above get error "Error: geom_line requires
    ## the following missing aesthetics: y" when the plot is finally
    ## printed at the end. Who knows why... even though printing the
    ## data frame and the computations below are apparently not
    ## affected by this at all (does not happen with the boxplot for
    ## simulated degree distribution)
    obs_deg_df$degree <- as.factor(obs_deg_df$degree)
    obs_deg_df$count[which(is.na(obs_deg_df$count))] <- 0
    obs_deg_df$nodefraction <- obs_deg_df$count / vcount(g_obs)
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
    if (is.bipartite(g_obs)) {
      degreetype <- ifelse(btype, 'mode B degree', 'mode A degree')
    } else {
      degreetype <- 'degree'
    }
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
    p <- p + guides(x = guide_axis(check.overlap = TRUE))
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
##    btype:      igraph bipartite node type FALSE or TRUE, or NULL
##                if not bipartite. (Default NULL)
##
## Return value:
##    ggplot2 object to add to plot list
##
deg_hist_plot <- function(g_obs, sim_graphs, mode, use_log, btype=NULL) {
#    print('in deg_hist_plot...')#XXX seems to be only way to debug in R...
    start <- Sys.time()
    if (use_log) {
      if (is.bipartite(g_obs)) {
        dobs <- data.frame(degree = log(degree(g_obs, V(g_obs)[which(V(g_obs)$type == btype)], mode=mode)),
                           group = 'obs')
      } else {
        dobs <- data.frame(degree = log(degree(g_obs, mode=mode)),
                           group = 'obs')
      }        
    } else {
      if (is.bipartite(g_obs)) {
        dobs <- data.frame(degree = degree(g_obs, V(g_obs)[which(V(g_obs)$type == btype)], mode=mode),
                           group = 'obs')
      } else {
        dobs <- data.frame(degree = degree(g_obs, mode=mode),
                           group = 'obs')
      }          
    }
#    print(names(dobs))#XXX
#    print('about to get simdegrees...')#XXX seems to be only way to debug in R...
  ## get degrees of all simulated graphs in one histogram
  if (is.bipartite(g_obs)) {
    ## as.vector() and unlist() BOTH seems to be required, otherwise rbind() below crashes with error about wrong number of columns
    simdegrees <- as.vector(unlist(sapply(sim_graphs, function(g) degree(g, V(g)[which(V(g)$type == btype)], mode=mode))))
    } else {
      simdegrees <- as.vector(unlist(sapply(sim_graphs, function(g) degree(g, mode=mode))))
    }
    if (use_log) {
        dsim <- data.frame(degree = log(simdegrees), group = 'sim')
    } else {
        dsim <- data.frame(degree = simdegrees, group = 'sim')
    }
#    print(names(dsim))#XXX
#    print('about to rbind dobs and dsim...') #XXX seems to be only way to debug in R...
    dat <- rbind(dobs, dsim)
    end <- Sys.time()
    cat(mode, "-degree histogram data frame construction took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    start <- Sys.time()
    ## https://stackoverflow.com/questions/29287614/r-normalize-then-plot-two-histograms-together-in-r
    p <- ggplot(dat, aes(degree, fill = group, colour = group)) +
        geom_histogram(aes(y = ..density..),
                       alpha = 0.4, position = 'identity', lwd = 0.2)
    if (is.bipartite(g_obs)) {
      degreetype <- ifelse(btype, 'mode B degree', 'mode A degree')
    } else {
      degreetype <- 'degree'
    }
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



##
## Return list of goodness-of-fit plots
##
## Parameters:
##    g_obs:       observed graph igraph object
##    sim_graphs:  simulated graphs list of igraph objects
##    do_subplots: if TRUE, do subplots of triad census on separate .eps file 
##                 (for bipartite graph, cycle length distribution is done
##                 on seperate .eps file, if do_cycledist is TRUE).
##                 also (default FALSE)
##    do_geodesic: if TRUE, include geodesic distance distribution plot
##                 otherwise do not (default TRUE); useful as this can
##                 be unusuably slow especially on large networks
##    do_dsp :     if TRUE, include dyadwise shared partner distribution
##                 otherwise do not (default TRUE); useful as this can use
##                 huge amounts of memory and be unsuably slow, even on
##                 bipartite networks (no ESP as always zero) that aren't even
##                 that big (few hundred thousand nodes)
##    do_bipartite_cc : if TRUE and graph is bipartite, then do bipartite
##                      clustering coefficients using tnet package.
##                      Default FALSE (as is very slow)
##    do_cycledist : If TRUE, do cycle length distribution (up to
##                   MAX_CYCLEDIST). Default FALSE (as is slow).
##                   For bipartite, this replaces four-cycle distribution
##                   (since that it included within this).
##
## Return value:
##    list of ggplot2 objects
##
##
build_sim_fit_plots <- function(g_obs, sim_graphs, do_subplots=FALSE,
                                do_geodesic=TRUE, do_dsp=TRUE,
                                do_bipartite_cc=FALSE, do_cycledist=FALSE) {

  num_sim <- length(sim_graphs)
  plotlist <- list()

  cat('obs num nodes: ', vcount(g_obs), '\n')
  cat('sim num nodes: ', sapply(sim_graphs, vcount), '\n')

  

  if (is.directed(g_obs)) {
    ##
    ## In degree
    ##

    system.time(plotlist <- c(plotlist,
                              list(deg_distr_plot(g_obs, sim_graphs, 'in'))))

    print('about to do deg_hist_plot...')
    system.time(plotlist <- c(plotlist,
                              list(deg_hist_plot(g_obs, sim_graphs, 'in', FALSE))))

    system.time(plotlist <- c(plotlist,
                              list(deg_hist_plot(g_obs, sim_graphs, 'in', TRUE))))


    ##
    ## Out degree
    ##

    system.time(plotlist <- c(plotlist,
                              list(deg_distr_plot(g_obs, sim_graphs, 'out'))))

    system.time(plotlist <- c(plotlist,
                              list(deg_hist_plot(g_obs, sim_graphs, 'out', FALSE))))

    system.time(plotlist <- c(plotlist,
                              list(deg_hist_plot(g_obs, sim_graphs, 'out', TRUE))))

  } else {
    ##
    ## Degree
    ##

    if (is.bipartite(g_obs)) {
      system.time(plotlist <- c(plotlist,
                                list(deg_distr_plot(g_obs, sim_graphs,
                                                    'all', FALSE))))
      system.time(plotlist <- c(plotlist,
                                list(deg_distr_plot(g_obs, sim_graphs,
                                                    'all', TRUE))))      
    } else {
      system.time(plotlist <- c(plotlist,
                                list(deg_distr_plot(g_obs, sim_graphs, 'all'))))
    }

    if (is.bipartite(g_obs)) {
      system.time(plotlist <- c(plotlist,
                                list(deg_hist_plot(g_obs, sim_graphs, 'all', FALSE, FALSE))))

      system.time(plotlist <- c(plotlist,
                                list(deg_hist_plot(g_obs, sim_graphs, 'all', TRUE, FALSE))))
      system.time(plotlist <- c(plotlist,
                                list(deg_hist_plot(g_obs, sim_graphs, 'all', FALSE, TRUE))))

      system.time(plotlist <- c(plotlist,
                                list(deg_hist_plot(g_obs, sim_graphs, 'all', TRUE, TRUE))))
    } else {
      system.time(plotlist <- c(plotlist,
                                list(deg_hist_plot(g_obs, sim_graphs, 'all', FALSE))))

      system.time(plotlist <- c(plotlist,
                                list(deg_hist_plot(g_obs, sim_graphs, 'all', TRUE))))
    }
  }

  ##
  ## Reciprocity
  ##
  if (is.directed(g_obs)) {
    ## Uses the default reciprocity in igraph which is probability
    ## that opposite counterpart of directed edge is also in the graph.
    ## Note can also use triad census 102 which is just graph with a mutual
    ## arc between two vertices (more related to the alternative reciprocity
    ## definition which we are not using), but not quite the same and also
    ## for very large graphs the triad 102 census count overflows and has to
    ## be omitted, while this does not.
    system.time( obs_reciprocity <- reciprocity(g_obs) )
    system.time( sim_reciprocity <- simplify2array(mclapply(sim_graphs, function(g) reciprocity(g))) )
    cat('obs reciprocity: ', obs_reciprocity, '\n')
    cat('sim reciprocity: ', sim_reciprocity, '\n')
    p <- ggplot() + geom_boxplot(aes(x = 'reciprocity', y = sim_reciprocity))
    p <- p + geom_point(aes(x = as.numeric(ordered('reciprocity')),
                            y = obs_reciprocity,
                            colour = obscolour))
    p <- p + ylab('fraction of arcs') + ptheme +
      theme(axis.title.x = element_blank())
    ##p <- p + ylim(0, 1)
    plotlist <- c(plotlist, list(p))
  }


  ##
  ## giant component size
  ##

  system.time(giant_component_sizes <- simplify2array(mclapply(sim_graphs,
                                              function(g) vcount(giant.component(g)))))
  giant_component_sizes <- giant_component_sizes / sapply(sim_graphs, vcount)
  obs_gcsize <- vcount(giant.component(g_obs)) / vcount(g_obs)
  cat('obs giant component size: ', obs_gcsize, '\n')
  cat('sim giant component size: ', giant_component_sizes, '\n')
  p <- ggplot() + geom_boxplot(aes(x = 'giant component', y = giant_component_sizes))
  p <- p + geom_point(aes(x = as.numeric(ordered('giant component')),
                          y = obs_gcsize,
                          colour = obscolour))
  p <- p + ylab('fraction of nodes')
  p <- p + ptheme +   theme(axis.title.x = element_blank())
  p <- p + theme(axis.text = element_text(colour = "black", size = rel(1.0))) # stop 'giant component' being light gray and smaller than other axis labels
  ##p <- p + ylim(0, 1)
  plotlist <- c(plotlist, list(p))



  ##
  ## Transitivity (global clustering coefficient and avg. local clustering coef.)
  ##
  
  cctypes <- c('average local', 'global') # must be in alphabetical order!
  system.time(ccs <- simplify2array(mclapply(sim_graphs, function(g) transitivity(g, type="global"))))
  system.time(cc_obs <- transitivity(g_obs, type='global'))
  cat('obs global cc: ', cc_obs, '\n')
  cat('sim global cc: ', ccs, '\n')
  system.time(ccs_localavg <- simplify2array(mclapply(sim_graphs, function(g)
    transitivity(g, type='localaverage'))))
  system.time(cc_localavg_obs <- transitivity(g_obs, type='localaverage'))
  cat('obs avg local cc: ', cc_localavg_obs, '\n')
  cat('sim avg local cc: ', ccs_localavg, '\n')
  if (is.bipartite(g_obs)) {
    cat('bipartite graph, not plotting clustering coefficients\n')
    # clustering coefficients must all be zero for bipartite graphs
    stopifnot(cc_obs == 0)
    stopifnot(all(ccs == 0))
    stopifnot(cc_localavg_obs == 0)
    stopifnot(all(ccs_localavg  == 0))
    
  } else {
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
    ##p <- p + ylim(0, 1)
    plotlist <- c(plotlist, list(p))
  }

  ##
  ## Triad census
  ##

  if (is.directed(g_obs)) {
    ## Note that on large networks, triad census counts can overflow and
    ## give negative numbers
    ## https://github.com/igraph/igraph/issues/625
    ## https://github.com/igraph/igraph/issues/497
    ## 003 and 012 and 102 can overflow on too large networks
    ## so will drop them if any are negative 
    dropFirstThree <- FALSE 
    nTriads_obs <- choose(vcount(g_obs), 3)
    nTriads_sim <- sapply(sim_graphs, function(h) choose(vcount(h), 3))
    system.time(obs_triadcensus <- triad.census(g_obs))
    num_triad_types <- length(obs_triadcensus)
    stopifnot(num_triad_types == 16)
    triadnames <- c('003', '012', '102', '021D', '021U', '021C', '111D',
                    '111U', '030T', '030C', '201', '120D', '120U', '120C',
                    '210', '300')
    stopifnot(length(triadnames) == num_triad_types)
    names(obs_triadcensus) <- triadnames
    cat('obs triad census: ', obs_triadcensus, '\n')
    if (obs_triadcensus[1] < 0 || obs_triadcensus[2] < 0 || obs_triadcensus[3] < 0) {
      dropFirstThree <- TRUE
    }
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
    obs_triadcensus_df$triadfraction <- obs_triadcensus_df$count / nTriads_obs
    for (i in 1:num_sim) {
      system.time(sim_triadcensus <- triad_census(sim_graphs[[i]]))
      names(sim_triadcensus) <- triadnames
      cat('sim triad census ', i, ': ', sim_triadcensus, '\n')
      if (sim_triadcensus[1] < 0 || sim_triadcensus[2] < 0 || sim_triadcensus[3] < 0) {
        dropFirstThree <- TRUE
      }
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

    sim_triadcensus_df$triadfraction <- sim_triadcensus_df$count / nTriads_sim
    ## Remove first two triads if necessary
    ## (003 triad (empty graph) and 012 (single edge))
    if (dropFirstThree) {
      cat ("WARNING: dropping triad census 003 and 012 and 102 as overflow detected\n")
      sim_triadcensus_df <- sim_triadcensus_df[which(sim_triadcensus_df$triad != "003"),]
      obs_triadcensus_df <- obs_triadcensus_df[which(obs_triadcensus_df$triad != "003"),]
      sim_triadcensus_df <- sim_triadcensus_df[which(sim_triadcensus_df$triad != "012"),]
      obs_triadcensus_df <- obs_triadcensus_df[which(obs_triadcensus_df$triad != "012"),]
      sim_triadcensus_df <- sim_triadcensus_df[which(sim_triadcensus_df$triad != "102"),]
      obs_triadcensus_df <- obs_triadcensus_df[which(obs_triadcensus_df$triad != "102"),]
    }
    p <- ggplot(sim_triadcensus_df, aes(x = triad, y = triadfraction))
    p <- p + geom_boxplot()
    p <- p + ylab('fraction of triads') + ptheme +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      xlab('Triad census')
    ## who knows why the hjust and vjust are needed, or what they should
    ## be, but they do seem to be, otherwise labels are not positioned right
    ## (note depends on which versoin of R/ggplot2 being used, but this worked
    ## when I wrote it with R 3.4.2 ggplot2 2.2.1 on Windows 10 cygwin:
    ## https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
    p <- p + geom_line(data = obs_triadcensus_df, aes(x = triad, y = triadfraction,
                                                      colour = obscolour,
                                                      group = 1))
    plotlist <- c(plotlist, list(p))  # no logarithm
    p <- p + scale_y_log10() + ylab("frac. triads (log scale)")
    plotlist <- c(plotlist, list(p))  # log scale on y axis

    if (do_subplots) {
      ## write separate file for triad census (log) plot
      ## add points for separate plot only (too large and messy on combined plots)
      ## also use raw number of triads rather than normalizing
      p <- ggplot(sim_triadcensus_df, aes(x = triad, y = count))
      p <- p + geom_boxplot()
      p <- p + ptheme +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        xlab('Triad census')
      ## who knows why the hjust and vjust are needed, or what they should
      ## be, but they do seem to be, otherwise labels are not positioned right
      ## (note depends on which versoin of R/ggplot2 being used, but this worked
      ## when I wrote it with R 3.4.2 ggplot2 2.2.1 on Windows 10 cygwin:
      ## https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
      p <- p + scale_y_log10(labels = my_scientific_10) 
      p <- p + ylab("number of triads (log scale)")
      p <- p + geom_line(data = obs_triadcensus_df, aes(x = triad, y = count,
                                                        colour = obscolour,
                                                        group = 1))
      p <- p + geom_point(data = obs_triadcensus_df, aes(x = triad, y = count,
                                                         colour = obscolour,
                                                         group = 1))
      triad_outfilename <- paste(simnetfileprefix, "_triadcensus.eps", sep="")
      cat("writing triad census (log) plot to EPS file ", triad_outfilename, "\n")
      postscript(triad_outfilename, horizontal=FALSE, onefile=FALSE, paper="special", width=9, height=6)
      ##pdf(triad_outfilename,onefile=FALSE, paper="special")
      print(p)
      dev.off()
    }


    ## ## log-odds version
    ## obs_triadcensus_df$logodds <- log(obs_triadcensus_df$triadfraction / (1 - obs_triadcensus_df$triadfraction))
    ## sim_triadcensus_df$logodds <- log(sim_triadcensus_df$triadfraction / (1 - sim_triadcensus_df$triadfraction))
    ## p <- ggplot(sim_triadcensus_df, aes(x = triad, y = logodds))
    ## p <- p + geom_boxplot()
    ## p <- p + ylab('log-odds') + ptheme +
    ##   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ##   xlab('Triad census')
    ## ## who knows why the hjust and vjust are needed, or what they should
    ## ## be, but they do seem to be, otherwise labels are not positioned right
    ## ## (note depends on which versoin of R/ggplot2 being used, but this worked
    ## ## when I wrote it with R 3.4.2 ggplot2 2.2.1 on Windows 10 cygwin:
    ## ## https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
    ## p <- p + geom_line(data = obs_triadcensus_df, aes(x = triad, y = logodds,
    ##                                                   colour = obscolour,
    ##                                                   group = 1))
    ## plotlist <- c(plotlist, list(p))  # no logarithm
  }


  ##
  ## geodesics (shortest paths)
  ##
  if (!do_geodesic) {
    cat("Not doing geodesic distance plot\n")
  } else if (vcount(g_obs) > MAX_SIZE_GEODESIC) {
    cat("WARNING: graph with ", vcount(g_obs),
        " too large to do geodesic fit, skipping\n")
  } else {
    num_dyads_obs <- choose(vcount(g_obs), 2) # N*(N-1)/2
    num_dyads_sim <- sapply(sim_graphs, function(h) choose(vcount(h), 2))
    system.time(obs_geodesics <- distance_table(g_obs)$res)
    system.time(sim_geodesics <- sapply(sim_graphs,
                                        function(g) distance_table(g)$res,
                                        simplify = FALSE))
    maxgeodesic <- max(length(obs_geodesics),
                       sapply(sim_geodesics, function(v) length(v)))
    cat("Max geodesic distance is ", maxgeodesic, "\n")
    geodesic_df <- data.frame(sim = rep(1:num_sim, each = maxgeodesic),
                              geodesic = rep(1:maxgeodesic, num_sim),
                              count = NA)
    start = Sys.time()
    for (i in 1:num_sim) {
      ## pad the sim vector to max length if it is not the longest already
      sg <- sim_geodesics[[i]]
      print(sg)#XXX
      print(length(sg))#XXX
      if (length(sg) < maxgeodesic) {
        oldlen <- length(sg)
        cat('oldlen = ', oldlen, '\n')#XXX
        sg <- rep(sg, length.out = maxgeodesic)
        sg[(oldlen+1):maxgeodesic] <- 0 # pad with zeroes
      }
      geodesic_df[which(geodesic_df[,"sim"] == i), "count"] <- sg
    }
    geodesic_df$geodesic <- as.factor(geodesic_df$geodesic)
    geodesic_df$nodefraction <- geodesic_df$count / num_dyads_sim
    end = Sys.time()
    cat("Geodesic sim data frame construction took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    start = Sys.time()
    ## pad the observed vector to max length if it is not the longest already
    if (length(obs_geodesics) < maxgeodesic) {
      oldlen <- length(obs_geodesics)
      obs_geodesics <- rep(obs_geodesics, length.out = maxgeodesic)
      obs_geodesics[(oldlen+1):maxgeodesic] <- 0 # pad with zeroes
    }
    obs_geodesic_df <- data.frame(geodesic = 1:maxgeodesic,
                                  count = as.numeric(obs_geodesics))
    obs_geodesic_df$nodefraction <- obs_geodesic_df$count / num_dyads_obs
    end = Sys.time()
    cat("Geodesic obs data frame construction took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    p <- ggplot(geodesic_df, aes(x = geodesic, y = nodefraction)) + geom_boxplot()
    p <- p + geom_line(data = obs_geodesic_df, aes(x = geodesic, y = nodefraction,
                                                   colour = obscolour, group = 1))
    p <- p + ptheme +
      xlab("geodesic distance") + ylab("fraction of dyads")
    p <- p + guides(x = guide_axis(check.overlap = TRUE))
    plotlist <- c(plotlist, list(p))
  }


  statnet_loaded <- FALSE

  ##
  ## edgewise shared partners
  ##

  if (is.bipartite(g_obs)) {
    cat("bipartite graph, not doing edgewise shared partners\n")
    # since ESP is always zero for bipartite graphs
  } else if (vcount(g_obs) > MAX_SIZE_ESP_DSP) {
    cat("WARNING: graph with ", vcount(g_obs),
        "nodes too large to do edgewise shared partners, skipping\n")
  } else {
    library(statnet)   
    library(intergraph)
    
    system.time(net_obs <- asNetwork(g_obs))
    system.time(sim_networks <- lapply(sim_graphs, function(g) asNetwork(g)))
    statnet_loaded <- TRUE # also net_obs and sim_networks built

    cutoff <- 50 # gw.cutoff default used in statnet is 30
    esp_df <- data.frame(sim = rep(1:num_sim, each = cutoff+1),
                         esp = rep(0:cutoff, num_sim),
                         count = NA)
    system.time(obs_esp <- summary(net_obs ~ esp(0:cutoff)))
    start <- Sys.time()
    for (i in 1:num_sim) {
      esp_df[which(esp_df[, "sim"] == i), "count"] <-  summary(sim_networks[[i]] ~ esp(0:cutoff))
      esp_df$edgefraction <- esp_df$count / network.edgecount(sim_networks[[i]])
    }
    end <- Sys.time()
    cat("esp sim data frame construction took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    obs_esp_df <- data.frame(esp = rep(0:cutoff),
                             count = summary(net_obs ~ esp(0:cutoff)))
    obs_esp_df$edgefraction <- obs_esp_df$count / network.edgecount(net_obs)
    end <- Sys.time()
    cat("esp obs data frame construction took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    ## remove zero counts from the end (use only up to max nonzero count)
    maxesp_sim <- max(esp_df[which(esp_df$count > 0),]$esp)
    maxesp_obs <- max(obs_esp_df[which(obs_esp_df$count > 0),]$esp)
    cat("Max obs esp is ", maxesp_obs, " and max sim esp is ", maxesp_sim, "\n")
    maxesp <- max(maxesp_sim, maxesp_obs)
    esp_df <- esp_df[which(esp_df$esp <= maxesp),]
    print("just before obs_esp_df subset") ## XXX trying to find line with error in R (why can't it just print line numbers like every other language?)
    obs_esp_df <- obs_esp_df[which(obs_esp_df$esp <= maxesp),]
    print("after obs_esp_df subset") ## XXX trying to find line with error in R (why can't it just print line numbers like every other language?)    
    obs_esp_df$esp <- as.factor(obs_esp_df$esp)
    print("just before factor(esp_df$esp)") ## XXX trying to find line with error in R (why can't it just print line numbers like every other language?)
    esp_df$esp <- as.factor(esp_df$esp)
    p <- ggplot(esp_df, aes(x = esp, y = edgefraction)) + geom_boxplot()
    p <- p + geom_line(data = obs_esp_df, aes(x = esp, y = edgefraction,
                                              colour = obscolour, group = 1))
    p <- p + ptheme + xlab("edgewise shared partners") +
      ylab("fraction of edges")
    p <- p + guides(x = guide_axis(check.overlap = TRUE))
    plotlist <- c(plotlist, list(p))

    ## add log scale version
    p <- p + scale_y_log10() + ylab("frac. edges (log scale)")
    plotlist <- c(plotlist, list(p))
  }


  ##
  ## dyadwise shared partners
  ##

  if (!do_dsp) {
    cat("Not doing dyadwise shared partners distribution\n")
  } else if (vcount(g_obs) > MAX_SIZE_ESP_DSP) {
    cat("WARNING: graph with ", vcount(g_obs),
        "nodes too large to do dyadwise shared partners, skipping\n")
  } else {
    library(statnet)   
    library(intergraph)

    system.time(net_obs <- asNetwork(g_obs))
    system.time(sim_networks <- lapply(sim_graphs, function(g) asNetwork(g)))
    statnet_loaded <- TRUE # also net_obs and sim_networks built

    cutoff <- 50 # gw.cutoff default used in statnet is 30
    dsp_df <- data.frame(sim = rep(1:num_sim, each = cutoff+1),
                         dsp = rep(0:cutoff, num_sim),
                         count = NA)
    system.time(obs_dsp <- summary(net_obs ~ dsp(0:cutoff)))
    start <- Sys.time()
    for (i in 1:num_sim) {
      dsp_df[which(dsp_df[, "sim"] == i), "count"] <-  summary(sim_networks[[i]] ~ dsp(0:cutoff))
      dsp_df$dyadfraction <- dsp_df$count / network.dyadcount(sim_networks[[i]])
    }
    end <- Sys.time()
    cat("dsp sim data frame construction took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    obs_dsp_df <- data.frame(dsp = rep(0:cutoff),
                             count = summary(net_obs ~ dsp(0:cutoff)))
    obs_dsp_df$dyadfraction <- obs_dsp_df$count / network.dyadcount(net_obs)
    end <- Sys.time()
    cat("dsp obs data frame construction took",
        as.numeric(difftime(end, start, unit="secs")), "s\n")
    ## remove zero counts from the end (use only up to max nonzero count)
    maxdsp_sim <- max(dsp_df[which(dsp_df$count > 0),]$dsp)
    maxdsp_obs <- max(obs_dsp_df[which(obs_dsp_df$count > 0),]$dsp)
    cat("Max obs dsp is ", maxdsp_obs, " and max sim dsp is ", maxdsp_sim, "\n")
    maxdsp <- max(maxdsp_sim, maxdsp_obs)
    dsp_df <- dsp_df[which(dsp_df$dsp <= maxdsp),]
    obs_dsp_df <- obs_dsp_df[which(obs_dsp_df$dsp <= maxdsp),]
    obs_dsp_df$dsp <- as.factor(obs_dsp_df$dsp)
    dsp_df$dsp <- as.factor(dsp_df$dsp)
    p <- ggplot(dsp_df, aes(x = dsp, y = dyadfraction)) + geom_boxplot()
    p <- p + geom_line(data = obs_dsp_df, aes(x = dsp, y = dyadfraction,
                                              colour = obscolour, group = 1))
    p <- p + ptheme + xlab("dyadwise shared partners") +
      ylab("fraction of dyads")
    p <- p + guides(x = guide_axis(check.overlap = TRUE))
    plotlist <- c(plotlist, list(p))
    ## add log scale version
    p <- p + scale_y_log10() + ylab("frac. dyads (log scale)")
    plotlist <- c(plotlist, list(p))
  }

  ##
  ## If the observed graph has self-edges (loops) then include those
  ##
  system.time( obs_loops <- sum(which_loop(g_obs)) )
  if (obs_loops > 0) {
    system.time( sim_loops <- sapply(sim_graphs, function(g) sum(which_loop(g))) )
    cat('obs loops: ', obs_loops, '\n')
    cat('sim loops: ', sim_loops, '\n')
    p <- ggplot() + geom_boxplot(aes(x = 'loops', y = sim_loops))
    p <- p + geom_point(aes(x = as.numeric(ordered('loops')),
                            y = obs_loops,
                            colour = obscolour))
    p <- p + ylab('count') + ptheme +
      theme(axis.title.x = element_blank())
    plotlist <- c(plotlist, list(p))
  }

  ##
  ## Four-cycles, for bipartite graphs only.
  ## Note that in a bipartite graph there are no 3-cycles, so 4-cycles
  ## are all chordless (induced) cycles.
  ##
  if (igraph::is.bipartite(g_obs) && !do_cycledist) {
    ## count_subgraph_isomorphisms() overcouts n-cycles by 2n as there are 2n
    ## automorphisms of an n-cycle.
    ## in igraph 1.3.0 could also use (faster) motifs(), but only have
    ## igraph 1.2.10 installed and not going through pain of 'updatings' just
    ## for this. See:
    ## https://stackoverflow.com/questions/71771349/counting-4-and-6-cycles-in-bipartite-r-igraph
    ## Note could also use statnet summary(g ~ cycle(4)) but that involves
    ## attaching statnet package which can cause huge memory probelems and
    ## also needing intergraph to convert (as required for ESP and DSP, see
    ## code above, which can be too slow and/or too much memory to work
    ## on large graphs).
    cat("bipartite graph, counting four-cycles...")
    ## Using igraph count_subgraph_isomorphisms() is too slow or fails
    ## with ' At vector.pmt:132 : cannot init vector, Out of memory'
    ## for networks only with a few thousand nodes, while statnet
    ## summary(g, ~ cycle(4)) only takes 5 seconds on the same network...
    ###system.time(obs_4cycles <- count_subgraph_isomorphisms(make_ring(4), g_obs)/
    ###                           (2*4))
    ###system.time(sim_4cycles <- sapply(sim_graphs, function(g)
    ###               count_subgraph_isomorphisms(make_ring(4), g) / (2*4)))
    if (!statnet_loaded) {
      library(statnet)
      library(intergraph)

      system.time(net_obs <- asNetwork(g_obs))
      system.time(sim_networks <- lapply(sim_graphs, function(g) asNetwork(g)))
      statnet_loaded <- TRUE # also net_obs and sim_networks built
    }

    print(system.time(obs_4cycles <- summary(net_obs ~ cycle(4))))
    print(system.time(sim_4cycles <- unlist(lapply(sim_networks,
                                      function(g) summary(g ~ cycle(4))))))
    cat('obs 4-cycles: ', obs_4cycles, '\n')
    cat('sim 4-cycles: ', sim_4cycles, '\n')
    p <- ggplot() + geom_boxplot(aes(x = 'four-cycles', y = sim_4cycles))
    p <- p + geom_point(aes(x = as.numeric(ordered('four-cycles')),
                            y = obs_4cycles,
                            colour = obscolour))
    p <- p + ylab('count') + ptheme +
      theme(axis.title.x = element_blank())
  p <- p + theme(axis.text = element_text(colour = "black", size = rel(1.0))) # stop 'four-cycles' being light gray and smaller than other axis labels
    plotlist <- c(plotlist, list(p))

  }

  ##
  ## Cycle length distribution (up to MAX_CYCLLEN only, as then gets 
  ## extremely slow and huge numbers of cycles).
  ##
  if (do_cycledist) {
    cat('cycle length distribution, MAX_CYCLELEN = ', MAX_CYCLELEN, '\n')
    if (!statnet_loaded) {
      library(statnet)
      library(intergraph)

      system.time(net_obs <- asNetwork(g_obs))
      system.time(sim_networks <- lapply(sim_graphs, function(g) asNetwork(g)))
      statnet_loaded <- TRUE # also net_obs and sim_networks built
    }
    if (igraph::is.bipartite(g_obs)) {
      cyclelens <- seq(4, MAX_CYCLELEN, 2) #only even length cycles in bipartite
    } else {
      cyclelens <- seq(3, MAX_CYCLELEN)
    }
    cat('computing cycle length distirbution in observed graph...')
    print(system.time(obs_cycledist <- summary(net_obs ~ cycle(cyclelens))))
    cat('obs_cycledist = ', obs_cycledist, '\n')
    cat('computing cycle length distirbution in simulated graphs...')
    print(system.time(sim_cycledist <- mclapply(sim_networks,
                                      function(g) summary(g ~ cycle(cyclelens)))))
    #print(sim_cycledist)
    cyclelen_df <- data.frame(sim = rep(1:num_sim, each = length(cyclelens)),
                               cyclelen = rep(cyclelens, num_sim),
                               count = NA)
    for (i in 1:num_sim) {
      cyclelen_df[which(cyclelen_df[,"sim"] == i), "count"] <- sim_cycledist[[i]]
    }
    obs_cyclelen_df <- data.frame(cyclelen = cyclelens,
                                  count = obs_cycledist)
    print(cyclelen_df)#XXX
    cyclelen_df$cyclelen <- factor(cyclelen_df$cyclelen)
    obs_cyclelen_df$cyclelen <- factor(obs_cyclelen_df$cyclelen)
    p <- ggplot(cyclelen_df, aes(x = cyclelen, y = count)) + geom_boxplot()
    p <- p + geom_line(data = obs_cyclelen_df, aes(x = cyclelen, y = count,
                                                   colour = obscolour, group = 1))
    p <- p + ptheme + xlab('cycle length') + ylab('count')
    p <- p + guides(x = guide_axis(check.overlap = TRUE))
    plotlist <- c(plotlist, list(p)) # no logarithm
    p <- p + scale_y_log10() + ylab("count (log scale)")
    plotlist <- c(plotlist, list(p))  # log scale on y axis
    if (do_subplots) {
      ## also write to separate file (add points for separate plot only,
      ## too large and messay on combined plots)
      p <- p + geom_point(data = obs_cyclelen_df, aes(x = cyclelen, y = count,
                                                     colour = obscolour, group = 1))
      p <- p + scale_y_log10(labels = my_scientific_10) 
      cycledist_outfilename <- paste(simnetfileprefix, "_cycledist.eps", sep="")
      cat("writing cycle length distribution plot to EPS file ", cycledist_outfilename, "\n")
      postscript(cycledist_outfilename, horizontal = FALSE, onefile = FALSE,
                 paper = "special", width = 9, height = 6)
      print(p)
      dev.off()  
    }
  }

  ##
  ## Bipartite clustering coefficients, for bipartite graphs only.
  ## This uses the tnet package, citation:
  ##
  ##    Opsahl, T., 2009. Structure and Evolution of Weighted Networks.
  ##    University of London (Queen Mary College), London, UK, pp. 104-122.
  ##    Available at http://toreopsahl.com/publications/thesis/;
  ##    http://toreopsahl.com/tnet/
  ##
  ## Note that this can be very slow for large networks, so may have to be
  ## disabled (option do_bipartite_cc=FALSE) for large networks.
  ## Citation for Robins-Alexander bipartite clustering coefficient is:
  ##
  ##   Robins, G., & Alexander, M. (2004). Small worlds among interlocking 
  ##   directors: Network structure and distance in bipartite graphs.
  ##   Computational & Mathematical Organization Theory, 10(1), 69-94.
  ##
  ## and for Opsahl bipartite clustering coefficient is:
  ##
  ##   Opsahl, T. (2013). Triadic closure in two-mode networks:
  ##   Redefining the global and local clustering coefficients.
  ##   Social networks, 35(2), 159-167.
  ##
  if (do_bipartite_cc && igraph::is.bipartite(g_obs)) {
    library(tnet)
    system.time(tn_obs <- as.tnet(get.edgelist(g_obs),
                                  type="binary two-mode tnet"))
    system.time(sim_tn <- mclapply(sim_graphs,
                                 function(g) as.tnet(get.edgelist(g),
                                                  type="binary two-mode tnet")))
    print(system.time(obs_robinsalexander_cc <- reinforcement_tm(tn_obs)))
    print(system.time(sim_robinsalexander_cc <- unlist(mclapply(sim_tn,
                                             function(g) reinforcement_tm(g)))))
    cat('obs Robins-Alexander cc: ', obs_robinsalexander_cc, '\n')
    cat('sim Robins-Alexander cc: ', sim_robinsalexander_cc, '\n')

    print(system.time(obs_opsahl_cc <- clustering_tm(tn_obs)))
    print(system.time(sim_opsahl_cc <- unlist(mclapply(sim_tn,
                                             function(g) clustering_tm(g)))))
    cat('obs Opsahl cc: ', obs_opsahl_cc, '\n')
    cat('sim Opsahl cc: ', sim_opsahl_cc, '\n')

    cctypes <- c('Opsahl', 'Robins-Alexander') # must be in alphabetical order!
    p <- ggplot() + geom_boxplot(aes(x = factor('Opsahl', levels=cctypes),
                                     y = sim_opsahl_cc))
    p <- p + geom_point(aes(x = as.numeric(factor('Opsahl', levels=cctypes)),
                            y = obs_opsahl_cc,
                            colour = obscolour))
    p <- p + geom_boxplot(aes(x = factor('Robins-Alexander', levels=cctypes),
                              y = sim_robinsalexander_cc))
    p <- p + geom_point(aes(x = as.numeric(factor('Robins-Alexander', levels=cctypes)),
                            y = obs_robinsalexander_cc,
                            colour = obscolour))
    p <- p + ylab('bipartite clustering coeff.') + ptheme +
      theme(axis.title.x = element_blank())
    ##p <- p + ylim(0, 1)
    plotlist <- c(plotlist, list(p))
  }



  return(plotlist)
}
