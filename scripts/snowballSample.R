#
# File:    snowballSample.R
# Author:  Alex Stivala
# Created: November 2013
#
# 
# Functions for
# snowball sampling in a (large) network, retaining zone information
# for each sampled node.
#

library(igraph)

#
# read_graph_file() - read adjacency matrix in paralle spnet (matrix) format
#                     or Pajek format
# 
# Parameters:
#     filename - filename of graph file to read from
#     directed - TRUE for directed graph
#
# Return value:
#     igraph graph object built from data in file
#
read_graph_file <- function(filename, directed) {
  alltext <- readLines(filename)
  if (any(grepl('*matrix', alltext, fixed=TRUE))) {
      return (read_graph_matrix_file(filename, directed))
  } else {
      pajek_text <- alltext
      if (any(grepl('***This graph contains:****', alltext, fixed=TRUE))) {
          # remove extra lines written by PNet
          pajek_text <- alltext[1:which(alltext == '***This graph contains:****')-1]
      }
      # remove all the vertex lines that don't seem to work with igraph pajek format
      firstline <- pajek_text[1]
      pajek_text <- pajek_text[(which(grepl('^[*]', pajek_text[2:length(pajek_text)]))+1):length(pajek_text)]
      tmpfilename <- tempfile()
      write(firstline, file=tmpfilename)
      write(pajek_text, file=tmpfilename, append=TRUE)
      g <- read.graph(tmpfilename, format="pajek")
      unlink(tmpfilename)
      stopifnot(is.directed(g) == directed)
      return(g)
  }
}

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

#
# giant_component() - return largest connected component of the graph
# 
# Paramters:
#    graph - igraph to get giant componetn of
#
# Return value:
#    largest connected component of graph
#
giant.component <- function(graph) {
  cl <- clusters(graph)
  return(induced.subgraph(graph, which(cl$membership == which.max(cl$csize))))
}


#
# snowball_sample() - snowball sampling
#
# Parameters:
#    g - graph to sample from( igraph)
#    num_waves - number of snowball waves (NB for cnostistney with SPNet
#                1 is subtraced so putting 3 here means there are really
#                only two waves [in normal usage] but there are 2 different
#                zones (0 for seeds, 1 wave first wave, 2 for 2nd wave)
#    seeds - vector of seeds (node ids) to start snowball sample from
#
# Return value:
#    graph (igraph) snowball sampled from g with each node having 
#    a zone attribute for snowball sampling zone(0=seed, 1=first wave, etc.)
#

# All sorts of fancy things can be done in R with igraph e.g.:
#
# induced.subgraph(h, Reduce(union, neighborhood(h, waves_num, sample.int(vcount(h), num_seeds) )))
# 
# or (to get the zone numbers (waves) instead of spedifying on input as above):
#
# lapply(sample.int(vcount(g), size=num_seeds), function(r) {zones <- graph.bfs(h, root=r, order=F, rank=F, father=F,pred=F,succ=F,dist=T,unreachable=F)$dist; list(nodes=which(zones<=num_waves), zones=zones[which(zones<=num_waves)])} )
#
# but in the end it is actually easier to just explictly write the code
# explicitly as done here to get the zones for all nodes without running
# into problems wih having to do union of node numbers and working out what
# zone corresponds to which node etc. all of which becomes a problem with
# "neat" versions like above in comments.
#
snowball_sample <- function(g, num_waves, seeds) {
  V(g)$zone <- NA
  V(g)[seeds]$zone <- 0
  nodes <- seeds
  newnodes <- nodes
  for (i in 1:num_waves) {
    newnodes <- Reduce(union, lapply(newnodes, function(v) 
                         Filter(function(x) !(x %in% nodes), neighbors(g, v))))
    if (!is.null(newnodes)) {
        V(g)[newnodes]$zone <- i
        nodes <- union(nodes, newnodes)
    }
  }
  return(induced.subgraph(g, nodes))
}

#
# snowball_sample_from_digraph() - snowball sampling from directed graph
#
# Parameters:
#    g - graph to sample from( igraph)
#    num_waves - number of snowball waves (NB for cnostistney with SPNet
#                1 is subtraced so putting 3 here means there are really
#                only two waves [in normal usage] but there are 2 different
#                zones (0 for seeds, 1 wave first wave, 2 for 2nd wave)
#    seeds - vector of seeds (node ids) to start snowball sample from
#
# Return value:
#    graph (igraph) snowball sampled from g with each node having 
#    a zone attribute for snowball sampling zone(0=seed, 1=first wave, etc.)
#
# This version does 'snowball sampling' from a directed graph by
# pretending that all the edges are actually unidrected (i.e. we 'follow'
# edges regardless of their direction - if there is a node a in our sample
# and directed link b->a then b will be sampled as we regrad it in the
# neighbourhood of a here (ignoring the edge direction). Hence this is
# is not really 'snowball sampling' as the graph we are sampling on
# is not the same as the actual graph - it is transformed by ignoring
# the edge dirctions (or equivalently thinking of them all as going 'out'
# from our seeds and samples). We might call it 'dirty snowball sampling',
# perhaps.
# All we have to do to implement this is use as.undirected(g) instead of just g
# in the neighbors() function.

snowball_sample_from_digraph <- function(g, num_waves, seeds) {
  V(g)$zone <- NA
  V(g)[seeds]$zone <- 0
  nodes <- seeds
  newnodes <- nodes
  for (i in 1:num_waves) {
    newnodes <- Reduce(union, lapply(newnodes, function(v) 
                         Filter(function(x) !(x %in% nodes),
                                neighbors(as.undirected(g), v))))
    if (!is.null(newnodes)) {
        V(g)[newnodes]$zone <- i
        nodes <- union(nodes, newnodes)
    }
  }
  return(induced.subgraph(g, nodes))
}

# 
# write_graph_file() - write adjacency matrix in parallel spnet format
#
# Parameters:
#     filename - filename to write to (warning: overwritten)
#     g - igrpah graph object
#     write_header - if TRUE write Pajek header lines
#
# Return value:
#    None
#
write_graph_file <- function(filename, g, write_header=TRUE) {
  f <- file(filename, open="wt")
  if (write_header) {
    cat('*vertices ', vcount(g), '\n', file=f)
    cat('*matrix\n', file=f)
  }
  write.table(get.adjacency(g, sparse=F), file=f, col.names=F, row.names=F)
  close(f)
}

#
# write_zone_file() - write zone file in parallel spnet (Pajek .clu) file format
#
# Parameters:
#    filename - filename to write to (warning: overwritten)
#    zones  - vector of snowball zones (waves) 0..n (0 for seed node)
#             elemetn i of the vector corresponds to node i of graph (1..N)
#
# Return value:
#    None.
#
write_zone_file <- function(filename, zones) {
  f <- file(filename, open="wt")
  cat('*vertices ', length(zones), '\n', file=f)
  write.table(zones, file=f, row.names=F, col.names=F)
  close(f)
}

# 
# write_subactors_file() - write subacctors file in parallel spnet format
#
# Parameters:
#     filename - filename to write to (warning: overwritten)
#     g - igrpah graph object
#
# Return value:
#    None
#
# See documetnation of this file in snowball.c readSubactorsFile()
# (format written by showActorsFile()).
#
write_subactors_file <- function(filename, g) {
  # TODO get attributes, currently just writes file with no attributes
  num_bin_attr = 0
  num_cont_attr = 0
  num_cat_attr = 0
  f <- file(filename, open="wt")
  cat("* Actors ", vcount(g), "\n", file=f)
  cat("* Number of Binary Attributes = ", num_bin_attr, "\n", file=f)
  cat("* Number of Continuous Attributes = ", num_cont_attr, "\n", file=f)
  cat("* Number of Categorical Attributes = ", num_cat_attr, "\n", file=f)
  cat("Binary Attributes:\n", file=f)
  cat("id\n", file=f) # TODO binary attriubute names
  for (i in 1:vcount(g)) {
    cat(i, file=f)
    for (j in 1:num_bin_attr) {
      cat (" ", file=f)
      # TODO output bin attr j value for node i
    }
    cat("\n", file=f)
  }
  cat("Continuous Attributes:\n", file=f)
  cat("id\n", file=f) # TODO continuous attriubute names
  for (i in 1:vcount(g)) {
    cat(i, file=f)
    for (j in 1:num_cont_attr) {
      cat (" ", file=f)
      # TODO output cont attr j value for node i
    }
    cat("\n", file=f)
  }
  cat("Categorical Attributes:\n", file=f)
  cat("id\n", file=f) # TODO categorical attriubute names
  for (i in 1:vcount(g)) {
    cat(i, file=f)
    for (j in 1:num_cat_attr) {
      cat (" ", file=f)
      # TODO output cat attr j value for node i
    }
    cat("\n", file=f)
  }
  close(f)
}

