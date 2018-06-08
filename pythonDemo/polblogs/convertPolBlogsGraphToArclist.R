#!/usr/bin/Rscript
#
# File:    convertPolBlogsGraphToArclist.R
# Author:  Alex Stivala
# Created: December 2014
#
# Convert political blogs graph GML from Newman's website http://www-personal.umich.edu/~mejn/netdata/polblogs.zip
# L. A. Adamic and N. Glance, "The political blogosphere and the 2004 US Election", in Proceedings of the WWW-2005 Workshop on the Weblogging Ecosystem (2005).
# to adjacency matrix for use in PNet
#
# Writes to the files polblogs_arclist.txt and polblogs_catattr.txt
#

library(igraph)

g <- read.graph('polblogs.gml', format='gml')
summary( g)
g <- simplify(g, remove.loops=TRUE, remove.multiple=TRUE)
summary( g)
write.graph(g, file="polblogs_arclist.txt", format=c("pajek"))
attr_file <- file('polblogs_catattr.txt', open='wt')
cat('value\n', file=attr_file)
cat(V(g)$value, file=attr_file, sep='\n')
close(attr_file)

