#!/usr/bin/Rscript
#
# File:    convertMatrixToArclist.R
# Author:  Alex Stivala
# Created: November 2016
#
#
#
#
# Read PNet adjacency matrix network and and convert to format
# for EstimNetDirected estimation. Input is stdin, output
# to stdout.
#
# Usage:
# 
# Rscript convertMatrixToArclist.R 
#

library(igraph)

#
# main
#


adjmatrix <- as.matrix(read.table(file("stdin")))
g <- graph.adjacency(adjmatrix, mode="directed")
# for some reason need binary connection for write.graph
# https://stackoverflow.com/questions/7422575/how-to-write-raw-type-bytes-to-stdout
con <- pipe("cat", "wb")
write.graph(g, con, format="pajek")
flush(con)
close(con)

