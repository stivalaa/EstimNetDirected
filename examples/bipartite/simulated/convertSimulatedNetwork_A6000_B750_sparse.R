#!/usr/bin/env Rscript
#
# Read bipartite network simulated by BPNet and write in simple bipartite
# Pajek .net format for EStimNetDirected
#
source('../../../scripts/snowballSample.R')

g <- read_graph_file('sample-statistics-A6000_B750_sparse_sim100000000.txt')
g <- remove.edge.attribute(g, 'weight')
write.graph(g, 'bpnet_A6000_B750_sparse_sim100000000.net', format='pajek')
