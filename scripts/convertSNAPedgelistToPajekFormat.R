#!/usr/bin/Rscript
#
# File:    convertSNAPedgelistToPajekFormat.R
# Author:  Alex Stivala
# Created: October 2018
#
# Read a gzipped edge list 
# from SNAP https://snap.stanford.edu/data/
# and convert to Pajek format
#
# Reference for SNAP collection of data sets:
# 
#@misc{snapnets,
#  author       = {Jure Leskovec and Andrej Krevl},
#  title        = {{SNAP Datasets}: {Stanford} Large Network Dataset Collection},
#  howpublished = {\url{http://snap.stanford.edu/data}},
#  month        = jun,
#  year         = 2014
#}
#
# Usage:
# 
# Rscript convert SNAPedgelistToPajekFormat.R infilename.txt.gz
#
# Output files (WARNING overwritten):
#     infilename.txt
#
#

library(igraph)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  cat("Usage: convertSNAPedgelistToPajekFormat.R\n")
  quit(save="no")
}
infile <- args[1]
basefilename <- sub("(.+)[.].+", "\\1", basename(infile))
#outfilename <- paste(basefilename, 'txt', sep='.')
outfilename <- basefilename

edgelist <- read.table(gzfile(infile))
# have to add 1 as igraph cannot handle 0 as vertex id apparently
g <- graph.edgelist(as.matrix(edgelist)+1, directed=TRUE)

# label nodes with 'name' which has special meaning in igraph allows lookup
# just use node id (starting at 1) as name
V(g)$name <- 1:vcount(g)

summary(g)
# remove multiple and self edges
g <- simplify(g , remove.multiple = TRUE, remove.loops = TRUE)
summary(g)


write.graph(g, outfilename, format="pajek")

