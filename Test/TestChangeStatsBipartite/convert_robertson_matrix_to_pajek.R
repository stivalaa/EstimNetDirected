#!/usr/bin/env Rscript
##############################################################################
#
# File:    convert_robertson_matrix_to_pajek.R
# Author:  Alex Stivala
# Created: April 2022
#
# Convert # Robertson (1929) pollinator network (two-mode) downloaded 
# from the Interaction Web Database 
# (IWDB) [http://www.ecologia.ib.usp.br/iwdb/index.html] 19 April 2022.
#  http://www.ecologia.ib.usp.br/iwdb/html/robertson_1929.html
# from biadjacency matrix format (rows pollinators, columns plants)
# to Pajek bipartite edgelist format
#
# Usage: Rscript convert_robertson_matrix_to_pajek.R
# 
# Reads data from data/ subdirector
# Writes output to ./robertson_pollinators_bipartite.net (WARNING: overwrites)
#
##############################################################################
library(igraph )

B <- as.matrix(read.table('data/robertson_1929_matr.txt', header=F))
m = dim(B)[1]
n = dim(B)[2]

cat('m = ', m, '\n')
cat('n = ', n, '\n')

# Convert biadjacency matrix to adjacency matrix
A <- rbind(cbind(matrix(0, m, m), B), cbind(t(B), matrix(0, n,n)))
g <- graph_from_adjacency_matrix(A, mode='undirected')
V(g)[1:m]$type <- FALSE
V(g)[(m+1):(m+n)]$type <- TRUE

summary(g)
summary(degree(g, V(g)[which(V(g)$type == FALSE)]))
summary(degree(g, V(g)[which(V(g)$type == TRUE)]))

write.graph(g, 'robertson_pollinators_bipartite.net', format='pajek')

