##
## File:    count_unique_bipartite_nodes_in_fourcycles.R
## Author:  Alex Stivala
## Created: October 2024
##
## Function to get lists of nodes in four-cycles in a graph.
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
## Usage: Rscript count_unique_bipartite_nodes_in_fourcycles.R <pajekfile>
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
  cat("Usage: Rscript count_unique_bipartite_nodes_in_fourcycles.R pajekfile.net\n", file=stderr())
  quit(save="no", status=1)
}
netfilename <- args[1]

g <- read.graph(netfilename, format='pajek')
print(summary(g))
if (!is.bipartite(g)) {
  cat('ERROR: graph is not bipartite\n')
  quit(save='no', status=1)
}

fourcycles_list <- get_fourcycle_nodes(g)
cat('Number of four-cycles:', length(fourcycles_list),'\n')
fourcycle_nodes <- unlist(fourcycles_list)
stopifnot(length(fourcycle_nodes) %% 4 == 0)
## must be same number of type FALSE(0) nodes as type TRUE(1) nodes in 4cycles
stopifnot(sum((V(g)[fourcycle_nodes]$type)) == length(fourcycle_nodes)/2)

unique_fourcycle_nodes <- unique(fourcycle_nodes)

## FALSE (0) type is A, TRUE (1) type is B
cat('Unique A nodes in four-cycles:', sum(1-V(g)[unique_fourcycle_nodes]$type), '\n')
cat('Unique B nodes in four-cycles:', sum(V(g)[unique_fourcycle_nodes]$type), '\n')


