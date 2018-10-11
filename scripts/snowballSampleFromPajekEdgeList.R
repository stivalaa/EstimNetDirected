#!/usr/bin/Rscript
##
## File:    snowballSampleFromEdgeList.R
## Author:  Alex Stivala
## Created: November 2013
##
## 
## Do snowball sampling in a (large) network, retaining zone information
## for each sampled node.
##
## Input file is Pajek format edge list.
##
## The graph may be directed or undirected. If directed we do snowball
## sampling on the unidrected version of the graph (i.e. ignore edge
## directions), and the sampled graph is the directed subgraph of the
## original directed graph induced by the nodes thus sampled.
##
## Output files (sample description file giving names of following
## files, subgraphs in Pajek format, zone files in Pajek .clu format
## giving zone for each node, attirbute files giving attributes for
## each node) in a directory in format used by EstimNetDirected
##
## Usage:
## 
## Rscript snowballSampleFromEdgelist.R [-d] adjlist.net num_samples num_seeds num_waves outputdir
##
##    -d : graph is directed, otherwise undirected
##    adjlist.net is Pajek format edge list as above
##    num_samples is number of snowball samples to create
##    num_seeds it number of seeds in each sample
##    num_Waves is numer of snowball sampling waves
##    outputdir is output directory to create output files in
##

library(igraph)

## read in R source file from directory where this script is located
##http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep=.Platform$file.sep))
}

source_local('snowballSample.R')


## 
## main
##


args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5 || length(args) > 6) {
  cat("Usage: snowballSample  [-d] edgelist.csv num_samples num_seeds num_waves  outputdirname\n")
  quit(save="no")
}
basearg <- 0
directed <- FALSE
if (length(args) == 6) {
    if (args[1] == "-d") {
        directed <- TRUE
        basearg <- 1
    }
    else {
        cat("Usage: snowballSample  [-d] edgelist.csv num_samples num_seeds num_waves  outputdirname\n")
        quit(save="no")
    }
}
input_pajek_file <- args[basearg+1]
num_samples <- as.integer(args[basearg+2])
num_seeds <- as.integer(args[basearg+3])
num_waves <- as.integer(args[basearg+4])
output_dir <- args[basearg+5]

if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

cat("directed: ", directed, "\n")
cat("number of samples: ", num_samples, "\n")
cat("number of seeds: ", num_seeds, "\n")
cat("number of waves: ", num_waves,"\n")
cat("output direcctory: ", output_dir, "\n")


cat("reading graph edge list from ", input_pajek_file, "...\n")
print(system.time(g <- read_graph_file(input_pajek_file, directed=directed)))
print(system.time(g <- simplify(g, remove.multiple=T, remove.loops=T)))
print(g) ## print info about graph

## get matrix where each row is random seed set for one sample
seedsets <- matrix(sample.int(vcount(g), num_seeds * num_samples, replace=F),
                   nrow = num_samples);

sampledesc_filename <- paste(output_dir, .Platform$file.sep,"sampledesc.txt",
                             sep="")
sampledesc_f <- file(sampledesc_filename, "wt")
for (i in 1:num_samples) {
  cat("generating snowball sample ", i, "...\n")
  if (directed) { 
      print(system.time(g_sample <- snowball_sample_from_digraph(g, num_waves,
                                                                 seedsets[i,])))
  }
  else {
      print(system.time(g_sample <- snowball_sample(g, num_waves, seedsets[i,])))
  }
  print(g_sample) ##print info about sample subgraph
  subgraph_filename <- paste(output_dir, .Platform$file.sep, "subgraph", i-1,
                             ".txt", sep="")
  subzone_filename <- paste(output_dir, .Platform$file.sep, "subzone", i-1,
                            ".clu", sep="")
  write_graph_file(subgraph_filename, g_sample)
  write_zone_file(subzone_filename,  V(g_sample)$zone)

  subactor_filename <- paste(output_dir, .Platform$file.sep, "subactor", i-1,
                             ".txt", sep="")
  ## TODO get actor attributes and write in EstimNetDirected format


  ## format of sampledesc file is:
  ## N subzone_filename subgraph_filename subactor_filename
  cat(vcount(g_sample), subzone_filename, subgraph_filename, subactor_filename,
             sep=" ", file=sampledesc_f)
  cat('\n', file=sampledesc_f)
}

close(sampledesc_f)

