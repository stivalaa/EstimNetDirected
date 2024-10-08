#!/usr/bin/Rscript
###
### File:    makeInouyePykeVisualizationAndStats.R
### Author:  Alex Stivala
### Created: August 2022
###
### Make visualization for Inouye-Pyke pollinator network
### and output some network statistics.
###
### Usage:
### 
### Rscript makeInouyePykeVisualizationAndStats.R 
###
### Output files (WARNING overwritten)
###    inouye_network.pdf      - visualization of biprtite network
###
### Network data citation:
###
###   Inouye, D. W., and G. H. Pyke. 1988. Pollination biology in the Snowy Mountains of Australia: 
###   comparisons with montane Colorado, USA. Australian Journal of Ecology 13:191-210.
###
### igraph citations:
###
###  Gabor Csardi (2015). igraphdata: A Collection of Network Data Sets
###  for the 'igraph' Package. R package version 1.0.1.
###  https://CRAN.R-project.org/package=igraphdata
###
###  Csardi G, Nepusz T: The igraph software package for complex network
###  research, InterJournal, Complex Systems 1695. 2006. http://igraph.org
###

library(igraph)


#
# read_graph_matrix_file() - read adjacency matrix in paralle spnet format
#
# Parameters:
#     filename - filename of subgraph file to read from
#     directed - TRUE for directed graph
#
# Return value:
#     igraph graph object built from adjacency matrix
#
read_graph_matrix_file <- function(filename, directed) {
  # skip all lines that do not start with 0 or 1
  # this includes the first two lines *vertices n and *matrix as well as
  # all lines following the actual matrix with various stats etc.
  # which SPnet writes there (but not snowballSample.R)
  alltext <- readLines(filename)
  matrixtext <- grep("^[01][01 ]+$", alltext, value=T)
  adjmatrix <- as.matrix(read.table(textConnection(matrixtext)))
  if (directed) {
      g <- graph.adjacency(adjmatrix, mode='directed')
  } else {
      g <- graph.adjacency(adjmatrix, mode='undirected')
  }
  return(g)
}


count_4cycles <- function(g) {
  ## https://lists.nongnu.org/archive/html/igraph-help/2009-04/msg00126.html
  ## https://stackoverflow.com/questions/71771349/counting-4-and-6-cycles-in-bipartite-r-igraph
  ## See also email "Gmail - biparite three-path change statistic"
  num_4cycles <- count_subgraph_isomorphisms(make_ring(4), g) / (2*4)
  return(num_4cycles)

###  ## in statnet this can be done much more easily and efficiently:
###  library(statnet)
###  library(intergraph)
###  gn <- asNetwork(g)
###  summary(gn ~ cycle(4))
}



############################################################################
##                          Main
############################################################################

set.seed(123)


args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 0) {
  cat("Usage: makeInouyePykeVisualizationAndStats.R\n")
  quit(save="no")
}


network_name <- 'inouye'

# read biadjacency matrix 
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

cat('four-cycles: ', count_4cycles(g), '\n')

print(g)

###
### write PDF file with network diagram
###
pdf('inouye_network.pdf')
plot(g, vertex.label = NA,
        vertex.size = 4,
        vertex.color = ifelse(V(g)$type, 'blue', 'red'),
        vertex.shape = ifelse(V(g)$type, 'square', 'circle'),
        layout = layout.auto
       )
dev.off()
