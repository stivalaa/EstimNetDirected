#!/usr/bin/env Rscript
##
## File:    convertStLouisCrimeDataToEstimNetDirectedFormat.R
## Auithor: Alex Stivala
## Created: April 2024
##
## Read St Louis crime network data from David Schoch's networkdata
## package (originally from Freeman data sets) and convert to
## EstimNetDirected format.
##
## Citation for data:
##
## S. Decker, C. W. Kohfeld, R. Rosenfeld, & J. Sprague, "St. Louis
## Homicide Project: Local Responses to a National Problem."
## University of Missouri-St. Louis (1991)
## https://books.google.com/books/about/The_St_Louis_Homicide_Project.html?id=umVAPQAACAAJ
##
## Note also available from archive of Freeman data sets (see below)
## in original UCINET format, and also from David Schoch's networkdata
## R package:
##
## https://github.com/schochastics/networkdata
##
## David Schoch. (2022). networkdata: Repository of Network Datasets
## (v0.1.14). Zenodo. https://doi.org/10.5281/zenodo.7189928
##
## and the networkdata version seems to line up exactly with the
## original data description (below), specifically by usnig the edge
## 'weight' as the type described below (1/2/3/4 =
## victim/suspect/witness/dual) rather than weight as such, while the
## KONECT version seems not:
##
## 	library(igraph)
## 	library(networkdata)
## 	data(crime)
## 	g<-upgrade_graph(crime)
##
## 	> g
## 	IGRAPH 58ccb2b UNWB 1427 1487 --
## 	+ attr: type (v/l), name (v/c), gender (v/n), weight (e/n)
## 	+ edges from 58ccb2b (vertex names):
## 	 [1] AbelDennis --950160 AbelDennis --960314 AbelDennis --960364
## 	 [4] AbelDennis --970473 AbramsChad --930029 AbramsChad --940121
## 	 [7] AbramsChad --940126 AbramsChad --940146 AbramsChad --950159
## 	[10] AbramsChad --950163 AbramsChad --950174 AbramsChad --950177
## 	[13] AbramsChad --950182 AbramsChad --950202 AbramsChad --950224
## 	[16] AbramsChad --960247 AbramsChad --960252 AbramsChad --960256
## 	[19] AbramsChad --960271 AbramsChad --960281 AbramsChad --960341
## 	[22] AbramsChad --960384 AbramsChad --960392 AbramsChad --970414
## 	+ ... omitted several edges
## 	>
## 	> sum(V(g)$type)
## 	[1] 557
## 	> summary(E(g)$weight)
## 	   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## 	  1.000   1.000   2.000   1.804   2.000   4.000
## 	> table(E(g)$weight)
##
## 	  1   2   3   4
## 	569 682 195  41
## 	>
##
##
## Original description from Freeman data sets moreno.ss.uci.edu,
## archived locally 2 October 2020 (original now unavailable) to
## file:///C:/Users/user/switchdrive/Institution/USI/datasets/FreemanData/moreno.ss.uci.edu/data.html#crime :
##
## =======================================================================
##
##  ROSENFELD,WHITE--ST. LOUIS CRIME
##
## DATASETS CRIME SEX
##
## DESCRIPTION One 870×557 two mode valued matrix of individuals by
## involvement in crime events. One 870×1 vector displaying the sex of
## each individual.
##
##     CRIME two mode, valued
##     SEX vector reporting sex of each individual. 
##
## BACKGROUND In the 1990s Rick Rosenfeld and Norm White used police
## records to collect data on crime in St. Louis. They began with five
## homicides and recorded the names of all the individuals who had
## been involved as victims, suspects or witnesses. They then explored
## the files and recorded all the other crimes in which those same
## individuals appeared. This snowball process was continued until
## they had data on 557 crime events. Those events involved 870
## participants of which: 569 appeared as victims 682 appeared as
## suspects 195 appeared as witnesses, and 41 were dual (they were
## recorded both as victims and suspects in the same crime. Their data
## appear, then, as an 870 by 557, individual by crime event
## matrix. Victims are coded as 1, suspects as 2, witnesses as 3 and
## duals as 4.  In addition Rosenfeld and White recorded the sex of
## each individual.
##
## =======================================================================
##
##
## Also add categorical attributes for male (1 for male, NA otherwise)
## and female (1 for female, NA otherwise).
## so that they can be used as categorical attributes with
## e.g. BipartiteNodematchBetaA in order to test differential homophily
## by having for each gender just one value and for the other NA so it
## never counts for matching (while the gender in question always matches
## itself only).
##
## Output files in cwd (WARNING: overwrites):
##     crime_bipartite.net
##     crime_binattr.net
##     crime_catattr.net
##
## Usage: Rscript convertStLouisCrimeDataToEstimNetDirectedFormat.R
##

library(igraph)
library(networkdata)
data(crime)

g <- crime
print(summary(g))
g <- remove.edge.attribute(g,"weight")
write.graph(g, 'crime_bipartite.net', format='pajek')
binattr <- data.frame(female =  ! V(g)$gender)
binattr$female <- ifelse(binattr$female, 1, 0)
print(summary(binattr))
write.table(binattr, 'crime_binattr.txt', row.names=F, col.names=T, quote=FALSE)
catattr <- data.frame(gender =  V(g)$gender) # female = 0, male=1
catattr$female <- ifelse(binattr$female == 1, 1, NA)
catattr$male <- ifelse(binattr$female == 0, 1, NA)
print(summary(catattr))
write.table(catattr, 'crime_catattr.txt', row.names=F, col.names=T, quote=FALSE)
