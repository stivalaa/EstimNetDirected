#!/usr/bin/env Rscript
##
## File:    count_fourcycles_per_bipartite_node_type.R
## Author:  Alex Stivala
## Created: April 2025
##
## Script to get histograms of counts of mode A and B nodes in four-cycles
##
## The CYPATH directory containing cypath executable and Perl scripts
## (transgrh.pl etc.) must be in the PATH as this R script uses
## system2() to call them in order to find (not just count)
## four-cycles. Have to edit the transgrh.pl script to add a
## semicolon on the end of line 38 as it does not work otherwise.
##
## For CYPATH see:
##
##    http://research.nii.ac.jp/~uno/code/cypath.html
##    http://research.nii.ac.jp/~uno/code/cypath11.zip
## 
##    Uno, T., & Satoh, H. (2014, October). An efficient algorithm for
##    enumerating chordless cycles and chordless paths. In International
##    Conference on Discovery Science (pp. 313-324). Springer, Cham.
##
##
## Usage: Rscript count_fourcycles_per_bipartite_node_type.R <pajekfile>
##
## Input file <pajekfie> is bipartite network in Pajek format.
## Output is to stdout.
##

library(igraph)

## read in R source file from directory where this script is located
##http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep=.Platform$file.sep))
}

source_local('get_fourcycle_nodes.R')



###
### Main
###

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  cat("Usage: Rscript count_fourcycles_per_bipartite_node_type.R pajekfile.net\n", file=stderr())
  quit(save="no", status=1)
}
netfilename <- args[1]

g <- read.graph(netfilename, format='pajek')
print(g)
if (!is.bipartite(g)) {
  cat('ERROR: graph is not bipartite\n')
  quit(save='no', status=1)
}
numnodesA <- sum(1-V(g)$type)
numnodesB <- sum(V(g)$type)
stopifnot(sum(1-V(g)$type) == vcount(g) - sum(V(g)$type))
cat('Number of A nodes: ', numnodesA, '\n')
cat('Number of B nodes: ', numnodesB, '\n')

fourcycles_list <- get_fourcycle_nodes(g)
cat('Number of four-cycles:', length(fourcycles_list),'\n')
fourcycle_nodes <- unlist(fourcycles_list)
stopifnot(length(fourcycle_nodes) %% 4 == 0)
## must be same number of type FALSE(0) nodes as type TRUE(1) nodes in 4cycles
stopifnot(sum((V(g)[fourcycle_nodes]$type)) == length(fourcycle_nodes)/2)

## FALSE (0) type is A, TRUE (1) type is B

c4distA <- sapply(which(!V(g)$type), function(i) sum(sapply(get_fourcycle_nodes(g), function(v) i %in% v)))
c4distB <- sapply(which(V(g)$type), function(i) sum(sapply(get_fourcycle_nodes(g), function(v) i %in% v)))
stopifnot(length(c4distA) == numnodesA)
stopifnot(length(c4distB) == numnodesB)
stopifnot(sum(c(c4distA, c4distB)) ==  4*length(fourcycles_list))

cat('Counts of four-cycles for mode A nodes:', c4distA, '\n')
cat('Counts of four-cycles for mode B nodes:', c4distB, '\n')



