#!/usr/bin/env Rscript
#
# Read bipartite network simulated by BPNet and write in simple bipartite
# Pajek .net format for EStimNetDirected
#
source('../../../scripts/snowballSample.R')

g <- read_graph_file('~/simulated_bipartite_networks/bpnet_A2000_B250_sim/sample-statistics-A2000_B250_sim5450000000.txt')
g <- remove.edge.attribute(g, 'weight')
write.graph(g, 'bpnet_A2000_B250_sim5450000000.net', format='pajek')
