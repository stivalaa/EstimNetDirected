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
MAX_SIZE_ESP_DSP <- 1000000 ## do not do shared partners if more nodes than this

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
##
## Return value:
##    ggplot2 object to add to plot list
##
##
deg_distr_plot <- function(g_obs, sim_graphs, mode) {
    num_sim <- length(sim_graphs)
    start = Sys.time()
    maxdeg <- max(unlist(sapply(sim_graphs, function(g) degree(g, mode=mode))),
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
    deg_df$nodefraction <- deg_df$count / sapply(sim_graphs, vcount)
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
##
## Return value:
##    ggplot2 object to add to plot list
##
deg_hist_plot <- function(g_obs, sim_graphs, mode, use_log) {
#    print('in deg_hist_plot...')#XXX seems to be only way to debug in R...
    start <- Sys.time()
    if (use_log) {
        dobs <- data.frame(degree = log(degree(g_obs, mode=mode)),
                           group = 'obs')
    } else {
        dobs <- data.frame(degree = degree(g_obs, mode=mode),
                           group = 'obs')
    }
#    print(names(dobs))#XXX
#    print('about to get simdegrees...')#XXX seems to be only way to debug in R...
    ## get degrees of all simulated graphs in one histogram
    simdegrees <- as.vector(unlist(sapply(sim_graphs, function(g) degree(g, mode=mode)))) ## as.vector() and unlist() BOTH seems to be required, otherwise rbind() below crashes with error about wrong number of columns
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



##
## Return list of goodness-of-fit plots
##
## Parameters:
##    g_obs:       observed graph igraph object
##    sim_graphs:  simulated graphs list of igraph objects
##    do_subplots: if TRUE, do subplots of triad census on separate .eps file 
##                 also (default FALSE)
##
## Return value:
##    list of ggplot2 objects
##
##
build_sim_fit_plots <- function(g_obs, sim_graphs, do_subplots=FALSE) {

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

    system.time(plotlist <- c(plotlist,
                              list(deg_distr_plot(g_obs, sim_graphs, 'all'))))

    system.time(plotlist <- c(plotlist,
                              list(deg_hist_plot(g_obs, sim_graphs, 'all', FALSE))))

    system.time(plotlist <- c(plotlist,
                              list(deg_hist_plot(g_obs, sim_graphs, 'all', TRUE))))
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
    system.time( sim_reciprocity <- sapply(sim_graphs, function(g) reciprocity(g)) )
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

  system.time(giant_component_sizes <- sapply(sim_graphs,
                                              function(g) vcount(giant.component(g))))
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
  ##p <- p + ylim(0, 1)
  plotlist <- c(plotlist, list(p))



  ##
  ## Transitivity (global clustering coefficient and avg. local clustering coef.)
  ##

  cctypes <- c('average local', 'global') # must be in alphabetical order!
  system.time(ccs <- sapply(sim_graphs, function(g) transitivity(g, type="global")))
  system.time(cc_obs <- transitivity(g_obs, type='global'))
  cat('obs global cc: ', cc_obs, '\n')
  cat('sim global cc: ', ccs, '\n')
  system.time(ccs_localavg <- sapply(sim_graphs, function(g)
    transitivity(g, type='localaverage')))
  system.time(cc_localavg_obs <- transitivity(g_obs, type='localaverage'))
  cat('obs avg local cc: ', cc_localavg_obs, '\n')
  cat('sim avg local cc: ', ccs_localavg, '\n')
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
  if (vcount(g_obs) > MAX_SIZE_GEODESIC) {
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
    plotlist <- c(plotlist, list(p))
  }



  ##
  ## edgewise shared partners
  ##

  if (vcount(g_obs) > MAX_SIZE_ESP_DSP) {
    cat("WARNING: graph with ", vcount(g_obs),
        "nodes too large to do edgewise shared partners, skipping\n")
  } else {
    library(statnet)   
    library(intergraph)
    
    system.time(net_obs <- asNetwork(g_obs))
    system.time(sim_networks <- lapply(sim_graphs, function(g) asNetwork(g)))

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

  if (vcount(g_obs) > MAX_SIZE_ESP_DSP) {
    cat("WARNING: graph with ", vcount(g_obs),
        "nodes too large to do dyadwise shared partners, skipping\n")
  } else {
    system.time(net_obs <- asNetwork(g_obs))
    system.time(sim_networks <- lapply(sim_graphs, function(g) asNetwork(g)))

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
  return(plotlist)
}
