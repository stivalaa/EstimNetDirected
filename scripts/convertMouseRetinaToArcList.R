#!/usr/bin/Rscript
#
# File:    convertMouseRetinaToArcList.R
# Author:  Alex Stivala
# Created: May 2019
#
# Read the GraphML file for mouse retina
# neurons downloaded from 
# https://neurodata.io/project/connectomes/
# and write Pajek arc list format.
#
# Citation for data is
#
#  M. Helmstaedter et al., "Connectomic reconstruction of the inner plexiform layer in the mouse retina." Nature 500, 168-174 (2013) Link
#
# Usage:
# 
# Rscript convertMouseRetinaToArcList.R 
#
# NB using simplify to remove multiedges and self-edges reduces
# original 577350 arcs to only 90811 arcs. This is entirely
# due to multiple edges; there are no self-loops.
#
# Output files (WARNING overwritten)
#    mouse_retina_arclist.txt         - arc list Pajek format
#    mouse_retina_catattr             - categorical attributes
#    mouse_retina_contattr            - continuous attributes
#
#

library(igraph)
# no longer seems to have if_else: library(dplyr) # massive overkill but needed for if_else since ifelse does not work for charcter data types

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 0) {
  cat("Usage: convertMouseRetinaToArcList.R\n")
  quit(save="no")
}

g <- read.graph("mouse_retina_1.graphml", format="graphml")
summary(g)
g <- simplify(g , remove.multiple = TRUE, remove.loops = TRUE)
summary(g)

##
## write graph
##
network_name <- "mouse_retina"
arclistfilename <- paste(network_name, "_arclist.txt", sep="")

write.graph(g, arclistfilename, format="pajek")
 
##
## write categorical attributes
##

catattr <- data.frame(designation = V(g)$designation,
                      volgyi_type = V(g)$volgyi_type,
                      macneil_type = V(g)$macneil_type)
# note have to use dplyr::if_else as ifelse forces everything into integer
# https://stackoverflow.com/questions/6668963/how-to-prevent-ifelse-from-turning-date-objects-into-numeric-objects
# why is everything in R always so difficult?
# also does not work as if_else in dplyr no longer seems to exist on the system I'm using??!?
# why is everying in R always so difficult?!?!!?!?
#catattr$designation <- if_else(catattr$designation == "NA", "NA", catattr$designation)
#catattr$volgyi_type <- if_else(catattr$volgyi_type == "NA" | catattr$volgyi_type == "", "NA", catattr$volgyi_type)
#catattr$macneil_type <- if_else(catattr$macneil_type == "NA" | catattr$macneil_type == "", "NA", catattr$macneil_type)
# so manually do it with lambda expression instead:
# (also note need to put as.character() everywhere, who knows why in R?)
setNA <- function(x) { if (x == "NA" | x == "") as.character(NA) else as.character(x) }
catattr$designation <- sapply(catattr$designation, setNA)
catattr$volgyi_type <- sapply(catattr$volgyi_type, setNA)
catattr$macneil_type <- sapply(catattr$macneil_type, setNA)

catattr$designation <- factor(catattr$designation)
print(levels(catattr$designation))
catattr$volgyi_type <- factor(catattr$volgyi_type)
print(levels(catattr$macneil_type))
catattr$macneil_type <- factor(catattr$macneil_type)
print(levels(catattr$macneil_type))
catattr$designation <- as.numeric(catattr$designation)
catattr$volgyi_type <- as.numeric(catattr$volgyi_type)
catattr$macneil_type <- as.numeric(catattr$macneil_type)
write.table(catattr, file = paste(network_name, "catattr.txt", sep="_"),
            row.names = FALSE, col.names = TRUE, quote=FALSE)

##
## write continuous attributes
##
contattr <- data.frame(x = V(g)$x,
                       y = V(g)$y,
                       z = V(g)$z)
write.table(contattr, file = paste(network_name, "contattr.txt", sep="_"),
            row.names = FALSE, col.names = TRUE, quote = FALSE)

