#!/usr/bin/env Rscript
##
## File:    convertOpsahlNewmanScientificCollaborationDataToEstimNetDirectedFormat.R
## Auithor: Alex Stivala
## Created: April 2024
##
## Read the Opsahl scientific collaboration data (Newman) bipartite
## network (author x article) and convert to format for EstimNetDirected.
##
## Citations for data:
##
## M. E. J. Newman, The structure of scientific collaboration networks,
## Proc. Natl. Acad. Sci. USA 98, 404-409 (2001).
##
## M. E. J. Newman, Scientific collaboration networks: I. Network
## construction and fundamental results, Phys. Rev. E 64, 016131 (2001).
##
## M. E. J. Newman, Scientific collaboration networks: II. Shortest paths,
## weighted networks, and centrality, Phys. Rev. E 64, 016132 (2001).
##
##  Opsahl, T. 2013. Triadic closure in two-mode networks: Redefining
##  the global and local clustering coefficients. Social Networks 35
##  (2), 159-167, doi: 10.1016/j.socnet.2011.07.001.
##
##
## Output files in cwd (WARNING: overwrites):
##     scientific_collaboration_bipartite.net
##
##
## Usage:
##    Rscript convertOpsahlNewmanScientificCollaborationDataToEstimNetDirectedFormat.R datadir
##    datadir is directory containing the downloaded data from
##    https://toreopsahl.com/datasets/#newman2001
##    Binary static two-mode network: tnet-format (659kb)
##    http://opsahl.co.uk/tnet/datasets/Newman-Cond_mat_95-99-two_mode.txt
##
##
## Information from https://toreopsahl.com/datasets/#newman2001 :
##
## 	Network 12: Newman’s scientific collaboration network
##
## 	This is the co-authorship network of based on preprints posted
## 	to Condensed Matter section of arXiv E-Print Archive between
## 	1995 and 1999. This dataset can be classified as a two-mode or
## 	affiliation network since there are two types of “nodes”
## 	(authors and papers) and connections exist only between
## 	different types of nodes. The two-mode network is projected
## 	onto one-mode networks using the procedure outlined on the
## 	projecting two-mode networks onto weighted one-mode
## 	networks-page. In addition to the network data, the names of
## 	the authors (369kb) are available.




library(igraph)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  cat("Usage: convertOpsahlNewmanScientificCollaborationDataToEstimNetDirectedFormat.R datadir\n")
  quit(save="no")
}
datadir <- args[1]


##
## read bipatrtite network in tnet format (two column edgelist; note that
## the identifiers for the two modes are not distinct)
##
dat <-read.table(file.path(datadir, 'Newman-Cond_mat_95-99-two_mode.txt'),
                 stringsAsFactors=FALSE, header=FALSE)

## Make string forms of identifiers so (a) the people and company identifiers
## are distinct, and (b) they can be used with igraph::make_bipartite_graph()
## aka igraph::graph.bipartite()

dat$AuthorID <- paste("Author", dat$V1, sep='.')
dat$ArticleID <- paste("Article", dat$V2, sep='.')

##
## Make logical node_types vector for make_bipartite_graph().
## It is a named vector where the names are the node names used
## in construting the bipartite graph.
## We will put all the persons (FALSE) first, followed by all the
## companies (TRUE).
## Note: needs igraph version higher than 1.2.11 to allow named nodes and
## types in make_bipartite_graph() (aka graph.bipartite()).
## Written using version 1.3.5.
## Bipartite network node attribute type is FALSE for persons (directors)
## and TRUE for companies (boards).
##
num_Authors <- length(unique(dat$AuthorID))
num_Articles <- length(unique(dat$ArticleID))
node_types <- c(rep(FALSE, num_Authors), rep(TRUE, num_Articles))
names(node_types) <- c(unique(dat$AuthorID), unique(dat$ArticleID))

##
## Build user x topic bipartite network
##
## Now make the edges vector in the correct format for make_bipartite_graph()
## I.e.:
##  "A vector defining the edges, the first edge points from the first
##   element to the second, the second edge from the third to the fourth, etc.
##   For a numeric vector, these are interpreted as internal vertex ids.
##   For character vectors, they are interpreted as vertex names. "
## [https://search.r-project.org/CRAN/refmans/igraph/html/make_graph.html]
##
## For the transpose and concatenate method to convert the two column
## edge list to this format, see
## https://stackoverflow.com/questions/41051823/convert-a-two-column-matrix-into-a-comma-separated-vector
##
gedges <- c(t(as.matrix(dat[,c("AuthorID", "ArticleID")])))
cat('num_Authors = ', num_Authors, '\n')
cat('num_Articles = ', num_Articles, '\n')
## make the bipartite graph
g <- graph.bipartite(types = node_types, edges = gedges, directed = FALSE)
print(summary(g))
stopifnot(typeof(V(g)$type) == 'logical')
stopifnot(sum(V(g)$type) == num_Articles)
stopifnot(num_Authors + num_Articles == vcount(g))
stopifnot(all(V(g)$type[1:num_Authors] == FALSE))
stopifnot(all(V(g)$type[(1+num_Authors):vcount(g)] == TRUE))


# there can be no self edges (since it is bipartite)
stopifnot(!any_loop(g))

## remove multiple edges (these can occur in this data as each row in the
## table is an appointment, so can be appointed to same board in difference
## capacities for example)
print('removing multiple edges...')
g <- simplify(g , remove.multiple = TRUE, remove.loops = FALSE)
summary(g)


##
## Writing graph in pajek format uses (as of when this 
## script was written, with igraph version 1.3.5) the internal
## integer node identifiers from 1..N which is just as 
## required for ALAAMEE or EstimNetDirected Pajek network format.
##
write.graph(g, 'scientific_collaboration_bipartite.net', format='pajek')


warnings()
