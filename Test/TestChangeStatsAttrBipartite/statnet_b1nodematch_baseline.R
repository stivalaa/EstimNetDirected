#!/usr/bin/Rscript
##
## File:    statnet_b1nodematch_baseline.R
## Author:  Alex Stivala
## Created: December 2023
##
## Get b1nodematch and b2nodematch statistics for test data using statnet
## as gold standard to verify implementations changeBipartiteNodematchBetaA,
## changeBipartiteNodematchBetaB, changeBipartiteNodematchAlphaA, and
## changeBipartiteNodematchAlphaB in EstimNetDirected.
##
library(igraph)
library(statnet)
library(intergraph)
g <- read.graph('../../examples/bipartite/simulated/bpnet_A12000_B4000_attrs_sim830000000.net', format='pajek')
catattr <- read.table('../../examples/bipartite/simulation/catattr_all.txt', header=T)
gn <- asNetwork(g)
gn <- as.network(as.edgelist(gn), bipartite=12000, directed=FALSE)
network::set.vertex.attribute(gn, 'catattrA', catattr$catattrA)
for (exponent in seq(0, 1, 0.1)) {
  print(summary(gn ~ b1nodematch("catattrA", alpha=exponent) + b1nodematch("catattrA", beta=exponent)))
}

