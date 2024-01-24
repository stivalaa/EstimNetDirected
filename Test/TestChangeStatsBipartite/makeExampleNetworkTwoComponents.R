#!/usr/bin/env Rscript
##
## File:    makeExampleNetworkTwoComponents.R
## Author:  Alex Stivala
## Created: November 2022
##
## Make Small biparite networks example with two components
## to demonstrate values of 
## four-cycle statistic vs BipartiteAltKCycleA and 
## BipartiteAltKCyclesB statistics.
##
## Usage: Rscript makeExampleNetworkTwoComponents.R
##
##  Writes .eps and .net files in cwd
##
##  WARNING: output files are overwritten
##
## Uses the igraph library to read Pajek format graph files and
## compute graph statistics:
##
##   Csardi G, Nepusz T: The igraph software package for complex network
##   research, InterJournal, Complex Systems
##   1695. 2006. http://igraph.org
## 
## Also uses the tnet package, citation:
##
##    Opsahl, T., 2009. Structure and Evolution of Weighted Networks.
##    University of London (Queen Mary College), London, UK, pp. 104-122.
##    Available at http://toreopsahl.com/publications/thesis/;
##    http://toreopsahl.com/tnet/
##
## Citation for Robins-Alexander bipartite clustering coefficient is:
##
##   Robins, G., & Alexander, M. (2004). Small worlds among interlocking 
##   directors: Network structure and distance in bipartite graphs.
##   Computational & Mathematical Organization Theory, 10(1), 69-94.
##
## and for Opsahl bipartite clustering coefficient is:
##
##   Opsahl, T. (2013). Triadic closure in two-mode networks:
##   Redefining the global and local clustering coefficients.
##   Social networks, 35(2), 159-167.
##

set.seed(12345)

library(igraph)
library(tnet)

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


plot_graph_viz <- function(g, outputepsfilename, layout=layout.auto) {
  postscript(outputepsfilename, onefile=FALSE)
  ##postscript(outputepsfilename, onefile=FALSE, fonts = "sans")
  plot(g,
  ##     vertex.label.family='sans',
       vertex.label = NA,
       vertex.color = ifelse(V(g)$type, 'blue', 'red'),
       vertex.shape = ifelse(V(g)$type, 'square', 'circle'),
       layout = layout
       )
  invisible(dev.off())
}

print_graph_stats <- function(g, name) {
## FIXME actually the incorrect values from reinforcement_tm() and
## clustering_tm() in tnet seem to be due to the conversion from
## igraph object to tnet object as.tnet(get.edgelist(g),...)
## here. It seems to work if intead the graph is written to Pajek
## file and then read back. Probably because tnet requires
## that the first column is mode A and the second mode B (resp. primary
## and secondary), while this is not what get.edgelist() in igraph returns
## (although it is this way if written and read back in Pajek format).
## So need to either do that or manipulate what get.edgelist() returns
## so that all the primary nodes are in first column and secondary nodes
## in second column.
## but for now have just manually computed Robins-Alexander bipartite c.c.
## as 4*C4/L3 using stats computed by EstimNetDirected - which in any
## case is far more efficient than tnet (but no stats for Opsahl bipartite
## c.c. implemented)

##  tn <- as.tnet(get.edgelist(g), type="binary two-mode tnet")
##  cat(name, vcount(g), ecount(g), count_4cycles(g),
##      reinforcement_tm(tn), clustering_tm(tn),
##     '\n')

  cat(name, vcount(g), ecount(g), count_4cycles(g), '\n')
}

##cat('Name', 'N', 'L', 'C4', 'Robins-Alexander', 'Opsahl', '\n')
cat('Name', 'N', 'L', 'C4', '\n')

g7 <- graph.ring(4) %>% add.vertices(1) %>% add.edges(c(1,5,3,5))
g9 <- g7 + graph.ring(4)
bimapping <- bipartite.mapping(g9)
stopifnot(bimapping$res)
V(g9)$type <- bimapping$type
print_graph_stats(g9, 'Four-cycles-4-two-components')
write.graph(g9,'fourcycle4_fourcycle_components_bipartite.net', format='pajek')
plot_graph_viz(g9, 'fourcycle4_fourcycle_components_bipartite.eps', layout = layout.sugiyama)
