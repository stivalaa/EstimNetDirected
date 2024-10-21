##
## File:    plot_fourcycles_dependence_figure.R
## Author:  Alex Stivala
## Created: June 2024
##
## Function to ake figure illustrating dependence structure for node-centered
## four-cycles statistic.
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


library(igraph)
library(RColorBrewer)

## read in R source file from directory where this script is located
##http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep=.Platform$file.sep))
}
source_local('get_fourcycle_edges.R')


##
## plot_fourcycles_dependence_figure() - plot figure showing 2-neighbourhoods
##
## 
## Parameters:
##     g           - igraph graph object to plot
##     i_node      - node in g for i of (i,j) dependency edge
##     j_node      - node in g for j of (i,j) dependency edge
##     outfilename - filename of PDF file to write (WARNING: overwrites)
##     k_node      - node to label as 'k', optional (default NULL)
##     l_node      - node to label as 'l', optional (default NULL)
##
## Return value:
##     None
##
plot_fourcycles_dependence_figure <- function(g, i_node, j_node, outfilename,
                                              k_node = NULL, l_node = NULL) {

  ## remove the i--j edge to compute neighbourhoods, it will be added back later
  g <- delete_edges(g, get.edge.ids(g, c(i_node, j_node)))
  
  bgroups <- c(ego(g, order=2, nodes=i_node),
               list(unique(unlist(ego(g, order=1,
                                      nodes = ego(g, order=1,
                                                  nodes=j_node)[[1]])))))
  names(bgroups) <- c("i", "j")

  ## https://www.r-bloggers.com/2018/05/visualizing-graphs-with-overlapping-node-groups/
  V(g)$color <- 'white'  
  group_colour <- brewer.pal(length(bgroups), 'Set1')
  group_colour_fill <- paste0(group_colour, '20')
  V(g)[bgroups[[1]]]$color <- group_colour_fill[1]
  V(g)[bgroups[[2]]]$color <- group_colour_fill[2]
  V(g)[i_node]$color <- group_colour[1]
  V(g)[j_node]$color <- group_colour[2]

  if (is.bipartite(g)) {
    V(g)$shape <- ifelse(V(g)$type, 'square', 'circle')
  }

  V(g)$label <- NA
  V(g)[i_node]$label <- "i"
  V(g)[j_node]$label <- "j"

  if (!is.null(k_node)) {
    V(g)[k_node]$label <- "k"
  }
  if (!is.null(l_node)) {
    V(g)[l_node]$label <- "l"
  }

  ## add back the i--j edge
  g <- g + edges(c(i_node, j_node))
  
  E(g)$lty<-1
  E(g)$color <- 'gray'
  E(g)[j_node %--% i_node]$lty <- 2

  ## nodevec <- unique(unlist(c(ego(g, order = 1, nodes = i_node),
  ##                            ego(g, order = 1, nodes = j_node))))
  ## nodevec <- unique(unlist(c(i_node,
  ##                            ego(g, order = 1, nodes = j_node))))
  if (is.bipartite(g)) {
    nodevec <- c(i_node)
  } else {
    nodevec <- c(i_node, j_node)
  }
  fourcycle_edges_list <- get_fourcycle_edges(g, nodes = nodevec)
  for (edgelist in fourcycle_edges_list) {
    E(g)[unlist(edgelist)]$color <- 'black'
  }

  print(outfilename)
  pdf(outfilename)
  par(mar = rep(0.1, 4))   # reduce margins
  plot(g, mark.groups=bgroups, mark.border=group_colour,
       mark.col = group_colour_fill, layout=layout.kamada.kawai, vertex.size=9,
       vertex.label.family = "sans", 
       vertex.label.font= 1) #plain text

  ## legend('topright',legend=names(bgroups),
  ##        col=group_colour,
  ##        pch = 18,
  ##        bty="n")
  dev.off()
}

