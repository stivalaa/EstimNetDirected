#!/usr/bin/env Rscript
##############################################################################
#
# File:    convert_inouye_pyke_matrix_to_pajek.R
# Author:  Alex Stivala
# Created: April 2022
#
# Convert Plant polinator web from Inouye & Pyke (1988) downloaded from IWDB 
# (Interaction Web DataBase).
# from biadjacency matrix format (rows pollinators, columns plants)
# to Pajek bipartite edgelist format

# Usage: Rscript convert_inouye_pyke_matrix_to_pajek.R
# 
# Reads data from data/ subdirector
# Writes output to ./inouye_pyke_pollinators_bipartite.net (WARNING: overwrites)
#
##############################################################################
library(igraph )

B <- as.matrix(read.table('data/inouye_matrix.txt', header=F))
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

write.graph(g, 'inouye_pyke_pollinators_bipartite.net', format='pajek')

