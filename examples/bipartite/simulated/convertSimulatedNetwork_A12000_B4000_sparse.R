#!/usr/bin/env Rscript
#
# Read bipartite network simulated by BPNet and write in simple bipartite
# Pajek .net format for EStimNetDirected
#
source('../../../scripts/snowballSample.R')

g <- read_graph_file('sample-statistics-A12000_B4000_sparse_sim770000000.txt')
g <- remove.edge.attribute(g, 'weight')
write.graph(g, 'bpnet_A12000_B4000_sparse_sim770000000.net', format='pajek')
