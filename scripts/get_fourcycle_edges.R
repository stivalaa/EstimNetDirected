##
## File:    get_fourcycle_edges.R
## Author:  Alex Stivala
## Created: June 2024
##
## Function to get lists of edges in four-cycles in a graph.
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

library(igraph)

##
## get_fourcycle_edges - return list of list of edges in 4-cycles
##
## Paramters:
##    g           - undirected graph object
##    nodes       - if not NULL, vector of nodes: only return cycles containing
##                  one or more of these nodes. Specified as node indices
##                  (integer starting at 1). Default NULL.
##
## Returns:
##    List of list where each vector is edges of g in cycle
##
## Note that the list of of lists (not vectors) of edges, by using
## SIMPLIFY = FALSE on mapply() [why is it SIMPLIFY on mapply()
## but simplify on sapply()? R is so confusing], because if
## SIMPLIFY = TRUE (the default) is used, then the igraph edge objects
## are simplified to just integers (index in E(g)) so less useful
## for e.g. checking if a particular edge specified by node pair
## is in it.
##
## Uses the CYPATH software to do this as there is no simple or
## efficient way in igraph.
##
## The directory of the CYPATH software from
## http://research.nii.ac.jp/~uno/code/cypath11.zip must be in
## thePATH. And the transgrh.pl script to add a semicolon on the end
## of lin#e 38 as it does not work otherwise.
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
get_fourcycle_edges <- function(g, nodes=NULL) {
  tmpfile <- tempfile()
  ## -1 to convert to zero-based for transgrh.pl and cypath
  eltext <- capture.output(write.table(get.edgelist(g, names=FALSE)-1,
                                       row.names = FALSE, col.names = FALSE))
  system2("transgrh.pl", input = eltext, stdout = tmpfile)

  ## Run cypath to list all four-cycles. The "c" option specifies
  ## cycles (could use "C" for chordless cycles - note in a bipartite
  ## graph all four-cycles are chordless).  The -l 4 and -u 4 options
  ## specify minimum and maximum, respectively, cycle length of 4.
  ## The "-, ," option sets the output delimiter to comma, so we can
  ## use grepl to get only lines with comma, so that summary (cycle
  ## length counts) information is removed (there is a q option
  ## documented to do this, but appears not to be implemented).
  
  cypath_output <- system2("cypath",
                           args = c("c", "-,", ",", "-l", "4", "-u", "4",
                                    tmpfile, "-"),
                           stdout = TRUE, stderr = FALSE)
  unlink(tmpfile)
  cycles <- cypath_output[grepl(',', cypath_output, fixed=TRUE)]
  ## +1 to convert back to one-based for igraph/R
  cycles <- lapply(cycles, function(s) as.numeric(strsplit(s, ',')[[1]]) + 1)

  if (!is.null(nodes)) {
    ## Only select cycles that contain a node in the nodes vector
    
    ## Note cypath has "-r n" to only include cycles contianing node n,
    ## but cannot do a list of nodes that way. TODO would be far more
    ## efficient if nodes is specified to cypath with -r multipel times
    ## rather than as currently where we get all cycles and just filter
    ## result.

    cycles <- Filter(function(x) any(nodes %in% x), cycles)
  }
  
  ## Now that each element of cycles list has vector of node
  ## identifiers in order along cycle, we can use consecutive nodes in
  ## each vector to get list of igraph edges.
  ## But first append first node on end (since it is a cycle) so we
  ## get every edge:
  cycles <- lapply(cycles, function(x) c(x, x[1]))
  cycle_edges <- lapply(cycles, function(x)
    mapply(function(i, j) E(g)[i %--% j], x[-length(x)], x[-1],
           SIMPLIFY = FALSE))
  
  return(cycle_edges)
}
