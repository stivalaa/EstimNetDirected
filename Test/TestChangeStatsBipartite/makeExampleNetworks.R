#!/usr/bin/env Rscript
##
## File:    makeExampleNetworks.R
## Author:  Alex Stivala
## Created: November 2022
##
## Make Small biparite networks examples to demonstrate values of 
## four-cycle statistic vs BipartiteAltKCycleA and 
## BipartiteAltKCyclesB statistics.
##
## Usage: Rscript makeExampleNetworks.R
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

g <- graph.ring(4) %>% add.vertices(1) %>% add.edges(c(1,5,3,5))
g <- g %>% add.vertices(1) %>% add.edges(c(1,6,3,6))
g <- g %>% add.vertices(1) %>% add.edges(c(1,7,3,7))
g <- g %>% add.vertices(1) %>% add.edges(c(1,8,3,8))
bimapping <- bipartite.mapping(g)
stopifnot(bimapping$res)
V(g)$type <- bimapping$type
print_graph_stats(g, 'Four-cycles-6')
write.graph(g,'fourcycle6_bipartite.net', format='pajek')
plot_graph_viz(g, 'fourcycle6_bipartite.eps', layout = layout.sugiyama)

g2 <- graph.ring(10)
stopifnot(bipartite.mapping(g2)$res)
V(g2)$type <- bipartite.mapping(g2)$type
print_graph_stats(g2, 'Ring')
write.graph(g2,'ring_bipartite.net', format='pajek')
plot_graph_viz(g2, 'ring_bipartite.eps')

g3 <- graph.lattice(dim = 1, length = 9)
stopifnot(bipartite.mapping(g3)$res)
V(g3)$type<-bipartite.mapping(g3)$type
print_graph_stats(g3, 'Chain')
write.graph(g3,'chain_bipartite.net', format='pajek')
plot_graph_viz(g3, 'chain_bipartite.eps')

g4 <- graph.star(10, mode="undirected")
stopifnot(bipartite.mapping(g4)$res)
V(g4)$type <- bipartite.mapping(g4)$type
print_graph_stats(g4, 'Star')
write.graph(g4,'star_bipartite.net', format='pajek')
plot_graph_viz(g4, 'star_bipartite.eps')

## Note this graph from Fig. 3 of Opsahl (2013) should have an Opsahl
## bipartite clustering coefficient of 0.6 (of 5 4-paths, 3 are closed;
## see p. 162 of opsahl (2013)). But we get NaN here, seems to be a bug?
## Even if it is due to 'primary' and 'secondary' node sets being
## different here, for a ring graph it should make no difference (symmetry)
## and yet we also get NaN for the ring graph. (And swappign the node
## types also made no difference; still get NaN)
## FIXME should we stop using the Opsahl bipartite c.c. as seems not
## correct in tnet implementation?
g5 <- graph.ring(6)
g5 <- g5 %>% add.vertices(1) %>% add.edges(c(2, 7))
stopifnot(bipartite.mapping(g5)$res)
V(g5)$type <- bipartite.mapping(g5)$type
##V(g5)$type <- !V(g5)$type # try swapping the node types
print_graph_stats(g5, 'Opsahl')
write.graph(g5,'opsahl_bipartite.net', format='pajek')
plot_graph_viz(g5, 'opsahl_bipartite.eps')

## It seems like reinforcement_tm() (Robins-Alexander bipartite c.c.)
## is also not correct: for a four-cycle get zero, but it should be 1.0
## as it is 4*C4/L3, and a four-cycle has C4=1 and L3=4
## 

g6 <- graph.ring(4)
bimapping <- bipartite.mapping(g6)
stopifnot(bimapping$res)
V(g6)$type <- bimapping$type
print_graph_stats(g6, 'Four-cycle')
write.graph(g6,'fourcycle_bipartite.net', format='pajek')
plot_graph_viz(g6, 'fourcycle_bipartite.eps', layout = layout.sugiyama)

g7 <- graph.ring(4) %>% add.vertices(1) %>% add.edges(c(1,5,3,5))
bimapping <- bipartite.mapping(g7)
stopifnot(bimapping$res)
V(g7)$type <- bimapping$type
print_graph_stats(g7, 'Four-cycles-3')
write.graph(g7,'fourcycle3_bipartite.net', format='pajek')
plot_graph_viz(g7, 'fourcycle3_bipartite.eps', layout = layout.sugiyama)

g8 <- graph.star(3, mode="undirected")
bimapping <- bipartite.mapping(g8)
stopifnot(bimapping$res)
V(g8)$type <- bimapping$type
print_graph_stats(g8, 'Two-path')
write.graph(g8,'twopath_bipartite.net', format='pajek')
plot_graph_viz(g8, 'twopath_bipartite.eps')

g10 <- graph.ring(6)
g10 <- g10 %>% add_edges(c(2,5))
bimapping <- bipartite.mapping(g10)
stopifnot(bimapping$res)
V(g10)$type <- bimapping$type
print_graph_stats(g10, 'Grid')
write.graph(g10,'grid_bipartite.net', format='pajek')
plot_graph_viz(g10, 'grid_bipartite.eps')

g11 <- g10 %>% delete.edges("5|6")
bimapping <- bipartite.mapping(g11)
stopifnot(bimapping$res)
V(g11)$type <- bimapping$type
print_graph_stats(g11, 'Grid-open')
write.graph(g11,'grid_open_bipartite.net', format='pajek')
plot_graph_viz(g11, 'grid_open_bipartite.eps')

g12 <- g7 %>% delete.edges("3|5")
bimapping <- bipartite.mapping(g12)
stopifnot(bimapping$res)
V(g12)$type <- bimapping$type
print_graph_stats(g12, 'Four-cycles-3-open')
write.graph(g12,'fourcycle3_open_bipartite.net', format='pajek')
plot_graph_viz(g12, 'fourcycle3_open_bipartite.eps', layout = layout.sugiyama)

g13 <- g7 %>% add_vertices(1) %>% add_edges(c(5,6))
bimapping <- bipartite.mapping(g13)
stopifnot(bimapping$res)
V(g13)$type <- bimapping$type
print_graph_stats(g13, 'Four-cycles-3-leaf')
write.graph(g13,'fourcycle3_leaf_bipartite.net', format='pajek')
plot_graph_viz(g13, 'fourcycle3_leaf_bipartite.eps', layout = layout.sugiyama)

g14 <- g13
V(g14)$type <- !V(g14)$type
print_graph_stats(g14, 'Four-cycles-3-leaf-revmode')
write.graph(g14,'fourcycle3_leaf_revmode_bipartite.net', format='pajek')
plot_graph_viz(g14, 'fourcycle3_leaf_revmode_bipartite.eps', layout = layout.sugiyama)

g15 <- g7
V(g15)$type <- !V(g15)$type
print_graph_stats(g15, 'Four-cycles-3-revmode')
write.graph(g15,'fourcycle3_revmode_bipartite.net', format='pajek')
plot_graph_viz(g15, 'fourcycle3_revmode_bipartite.eps', layout = layout.sugiyama)

##
## See email "New four-cycles ERGM parameter" sent 13 Feburary 2024:
##
##   So far I have only a single FourCyclesNodePower parameter (and is
##   equally applicable to one or two mode networks). But actually,
##   since it is based on per- node counts, it is easily separable
##   into FourCyclesNodePowerB and FourCyclesNodePowerA for the two
##   modes separately. Which for two-mode (bipartite) networks, is
##   likely to be desirable, e.g. for something like this structure
##   they have different values (only the central node's count is > 1
##   and hence affected by alpha) which is hidden if you just have the
##   single statistic ignoring mode. (As far as I know this is not a
##   named graph [https://en.wikipedia.org/wiki/List_of_graphs], but
##   it is something like a bipartite version of the "Dutch windmill
##   graph' otherwise known as the "Friendship graph" or "n-fan",
##   where instead of joining n copies of the triangle graph C3 to a
##   common vertex, I have joined n copies (n = 3 in the drawing) of
##   the four-cycle C4 (which in a bipartite graph is necessarily
##   chordless i.e. a hole) to a common vertex in a bipartite graph.)
##

g16 <- graph.ring(4) %>% add.vertices(3) %>%add.edges(c(1,5,5,6,6,7,7,1)) %>%
  add.vertices(3) %>%  add.edges(c(1,8,8,9,9,10,10,1))
bimapping <- bipartite.mapping(g16)
stopifnot(bimapping$res)
V(g16)$type <- bimapping$type
print_graph_stats(g16, 'Four-fan-3')
write.graph(g16,'fourfan3_bipartite.net', format='pajek')
plot_graph_viz(g16, 'fourfan3_bipartite.eps', layout = layout.kamada.kawai)


##
## Reconstructing examples from Peng (email reply to "4C -Paper" 31/1/2025)
##

## pw1 is actually the same as Four-cycles-3
pw1 <- graph.empty(n = 5, directed = FALSE)
V(pw1)$type <- c(rep(FALSE, 2), rep(TRUE, 3))
pw1 <- add.edges(pw1, c(1,3, 1,4, 1,5, 2,3, 2,4, 2,5))
stopifnot(is.bipartite(pw1))
print_graph_stats(pw1, 'pw1')
write.graph(pw1, 'pw1_bipartite.net', format='pajek')
plot_graph_viz(pw1, 'pw1_bipartite.eps', layout = layout.sugiyama)


pw2 <- graph.empty(n = 6, directed = FALSE)
V(pw2)$type <- c(rep(FALSE, 3), rep(TRUE, 3))
pw2 <- add.edges(pw2, c(1,4, 1,6, 2,4, 2,5, 2,6, 3,5, 3,6))
stopifnot(is.bipartite(pw2))
print_graph_stats(pw2, 'pw2')
write.graph(pw2, 'pw2_bipartite.net', format='pajek')
plot_graph_viz(pw2, 'pw2_bipartite.eps', layout = layout.sugiyama)


pw3 <- graph.empty(n = 7, directed = FALSE)
V(pw3)$type <- c(rep(FALSE, 4), rep(TRUE, 3))
pw3 <- add.edges(pw3, c(1,5, 1,7, 2,5, 2,7, 3,6, 3,7, 4,6, 4,7))
stopifnot(is.bipartite(pw3))
print_graph_stats(pw3, 'pw3')
write.graph(pw3, 'pw3_bipartite.net', format='pajek')
plot_graph_viz(pw3, 'pw3_bipartite.eps', layout = layout.sugiyama)
