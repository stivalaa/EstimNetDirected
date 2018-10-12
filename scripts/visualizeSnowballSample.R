#!/usr/bin/Rscript
##
## File:    visualizeSnowballSample.R
## Author:  Alex Stivala
## Created: November 2013
##
## Draw plot of snowball sample crated by snowballSampleFromNexus.R
## etc., read from the output format of those scripts (edge list in
## Pajek format format and zone attribute file for zone numbers). 
##
## Usage: Rscript visualizeSnowballSample.R [-d] subgraphfile.txt zonefile.clu output_postscript_file.eps
##
##  -d  : graph is directed (default undirected)
##  subgraphfile.txt is matrix file output of snowballSample.R write_graph_file()
##  zonefile.txt     is zone file output of snowballSample.R write_zone_file()
##  output_podf_file.eps is filename to write postscript of graph viz to (WARNING: overwritten)
##

library(igraph)
library(RColorBrewer)

## read in R source file from directory where this script is located
##http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep=.Platform$file.sep))
}

source_local('snowballSample.R')



## 
## read_zone_file() - read zone file in EstimNetDirected format
##                    It is just header line "zone" followed by zone number
##                    (0,1,2,...) of each node one per line in sequential
##                    nodeid order.
##
## Prameters:
##
##     filename - filename to read from
##
## Return value:
##     vector of zone numbers 
##
read_zone_file <- function(filename) {
   ## skip first line which is zone header line
   zone_df <- read.table(filename, header=F, skip=1)
   return(zone_df$V1)
}

## 
## main
##

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3 || length(args) > 4) {
  cat("Usage: Rscript visualizeSnowballSample.R [-d] graph_filename.txt zone_filename.txt output_file.eps\n")
  quit(save="no")
}
basearg <- 0
directed = FALSE
if (length(args) == 4) {
    if (args[1] == "-d") {
        directed <- TRUE
        basearg <- basearg + 1
    }
    else {
        cat("Usage: Rscript visualizeSnowballSample.R [-d] graph_filename.txt zone_filename.txt output_file.eps\n")
        quit(save="no")
    }
}
graph_filename <- args[basearg+1]
zone_filename <- args[basearg+2]
output_postscript_filename <- args[basearg+3]


g <- read_graph_file(graph_filename, directed)
zones <- read_zone_file(zone_filename)
V(g)$zone <- zones ## zone numbers line up with vertices
zones <- V(g)$zone
palette <- rev(brewer.pal(max(zones)+1, "YlGnBu")) ##colorblind-safe,print friendly
colours <- palette[(zones+1)]
## for inclusion in LaTeX or epstopdf, stop it rotating it
##postscript(output_postscript_filename,horizontal=FALSE,onefile=FALSE,paper='special',width=9,height=6)
postscript(output_postscript_filename,horizontal=FALSE,onefile=FALSE)
##http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
layout(rbind(1,2), heights=c(9,1)) ## put legend on bottom 1/10th of the chart
par(mar=c(0,0,0,0))
plot(g, layout=layout.fruchterman.reingold, vertex.size=5,
     ##vertex.label=V(g)$zone, vertex.label.family="Helvetica",
     vertex.label=NA,
     edge.arrow.size = 0.5,
     vertex.color = colours)

## setup for no margins on the legend
par(mar=c(0,0,0,0))
## c(bottom,left,top,right)
plot.new()
zonelabels <- sort(unique(zones))
legend("center", horiz=TRUE, legend=zonelabels, 
       pch=19, col=palette, title="Wave",
       ##pt.cex=1,
       cex=1.5)

dev.off()


