#!/usr/bin/env Rscript
#
# make viualizations of simulated graphs
#
# Uses igraph and graphlayouts
#
#  David Schoch (2020). graphlayouts: Additional Layout Algorithms for
#  Network Visualizations. R package version 0.7.1.
#  https://CRAN.R-project.org/package=graphlayouts
#
#
library(igraph)
library(graphlayouts)

##
## giant_component_and_rest() - return largest connected component of the graph
##                              as first element of list and the rest as
##                              second element of list
##
## Paramters:
##    graph - igraph to get giant componetn of
##
## Return value:
##    largest connected component of graph
##
giant_component_and_rest <- function(graph) {
  cl <- clusters(graph)
  return(list(
    gc = induced.subgraph(graph,
                          which(cl$membership == which.max(cl$csize))),
    rest = induced.subgraph(graph,
                            which(cl$membership != which.max(cl$csize)))
  ))

}



##
## Plotting the whole graph squashes giant component and takes a lot
## of space for the many smaller compoennts, which is messy and misleading
## since the giant component is 86% of the network. So instead plot
## the giant component and the rest of it separtely, but using par(mfrow(..))
## to put on smae plot.
##
##
plotgraph <- function(g) {
  plot(g, vertex.label=NA, vertex.size=2, vertex.frame.color = NA,
       vertex.shape=ifelse(V(g)$type,'square', 'circle'),
       vertex.color = V(g)$type+1,
       edge.width = 0.1,
       edge.color = "black",
       layout = layout_with_stress)
}


g <- read.graph('bpnet_A12000_B4000_sparse_sim770000000.net', format='pajek')
pdf('bpnet_A12000_B4000_sparse_sim770000000_viz.pdf')
par(mfrow=c(1, 2))
## reduce white space around igraph plot
## https://lists.nongnu.org/archive/html/igraph-help/2011-07/msg00036.html
par(mar=c(0,0,0,0)+.1)
comps <- giant_component_and_rest(g)
system.time(plotgraph(comps$gc))   # giant component
system.time(plotgraph(comps$rest)) # everything else
dev.off()

g <- read.graph('bpnet_A6000_B750_sparse_sim100000000.net', format='pajek')
pdf('bpnet_A6000_B750_sparse_sim100000000_viz.pdf')
par(mfrow=c(1, 2))
## reduce white space around igraph plot
## https://lists.nongnu.org/archive/html/igraph-help/2011-07/msg00036.html
par(mar=c(0,0,0,0)+.1)
comps <- giant_component_and_rest(g)
system.time(plotgraph(comps$gc))   # giant component
system.time(plotgraph(comps$rest)) # everything else
dev.off()

g <- read.graph('bpnet_A750_B250_sim777000000.net', format='pajek')
pdf('bpnet_A750_B250_sim777000000_viz.pdf')
par(mfrow=c(1, 2))
## reduce white space around igraph plot
## https://lists.nongnu.org/archive/html/igraph-help/2011-07/msg00036.html
par(mar=c(0,0,0,0)+.1)
comps <- giant_component_and_rest(g)
system.time(plotgraph(comps$gc))   # giant component
system.time(plotgraph(comps$rest)) # everything else
dev.off()
