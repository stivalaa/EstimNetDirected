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

source('../../../scripts/snowballSample.R')

g <- read.graph('bpnet_A12000_B4000_sparse_sim770000000.net', format='pajek')
pdf('bpnet_A12000_B4000_sparse_sim770000000_viz.pdf')
system.time( plot(g, vertex.label=NA,vertex.size=2, vertex.color = V(g)$type+1, layout=layout_with_stress) )
dev.off()
pdf('bpnet_A12000_B4000_sparse_sim770000000_gc_viz.pdf')
system.time( plot(giant.component(g), vertex.label=NA,vertex.size=2, vertex.color = V(giant.component(g))$type+1, layout=layout_with_stress, edge.width=.1,edge.color='black',vertex.frame.color=NA) )

g <- read.graph('bpnet_A6000_B750_sparse_sim100000000.net', format='pajek')
pdf('bpnet_A6000_B750_sparse_sim100000000_gc_viz.pdf')
system.time( plot(giant.component(g), vertex.label=NA,vertex.size=2, vertex.color = V(giant.component(g))$type+1, layout=layout_with_stress, edge.width=.1,edge.color='black',vertex.frame.color=NA) )
dev.off() 

g <- read.graph('bpnet_A750_B250_sim777000000.net', format='pajek')
pdf('bpnet_A750_B250_sim777000000_viz.pdf')
system.time( plot(g, vertex.label=NA,vertex.size=2, vertex.color = V(g)$type+1, layout=layout_with_stress, edge.width=.1,edge.color='black',vertex.frame.color=NA) )
dev.off()
pdf('bpnet_A750_B250_sim777000000_gc_viz.pdf')
system.time( plot(giant.component(g), vertex.label=NA,vertex.size=2, vertex.color = V(giant.component(g))$type+1, layout=layout_with_stress, edge.width=.1,edge.color='black',vertex.frame.color=NA) )
dev.off()
pdf('bpnet_A750_B250_sim777000000_gc_viz2.pdf')
system.time( plot(giant.component(g), vertex.label=NA,vertex.size=2, vertex.color = V(giant.component(g))$type+1, edge.width=.1,edge.color='black',vertex.frame.color=NA, layout = layout_with_centrality(giant.component(g), cent=betweenness(giant.component(g)))) )
dev.off()
