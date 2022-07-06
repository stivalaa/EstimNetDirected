##!/usr/bin/Rscript
##
## File:    convertEvtushenkoBoardDirectorAffiliationsToEstimNetDirectedFormat.R
## Author:  Alex Stivala
## Created: July 2022
##
## Read the network and attributes data in GML format for the 
## Evtushenko & Gastner (2019) company director network 
## and convert to Pajek bipartite format for EstimNetDirected. 
##
## 
## Usage: Rscript convertEvtushenkoBoardDirectorAffiliationsToEstimNetDirectedFormat.R  board_director.gml.gz
## where board_director.gml.gz is the file downloaded from 
## https://zenodo.org/record/3553442
## (publication date November 26, 2019; downloaded 28 June 2020)
#
## Output files in cwd (WARNING overwritten):
##     evtushenko_directors_bipartite.net
##     evtushenko_directors_binattr.txt
##     evtushenko_directors_catattr.txt
##     evtushenko_directors_contattr.txt
##
##
## Citation for publication of data:
##
##  A. Evtushenko and M. T. Gastner, Beyond Fortune 500: Women in a global
##  network of directors. In H. Cherifi et al. (Eds.), Complex Networks and 
##  Their Applications VIII, Proc. 8th Int. Conf. Complex Networks and Their 
##  Applications, Volume 1, pp. 586-598 (Springer, Cham, 2020), 
##  DOI: 10.1007/978-3-030-36683-4_47.
##
## Citation for data:
## 
##  Evtushenko, Anna, & Gastner, Michael T. (2019). Data set discussed in 
##  "Beyond Fortune 500: Women in a Global Network of Directors" (1.0)
##  [Data set]. 8th Int. Conf. Complex Networks and Their Applications
##  (COMPLEX NETWORKS 2019), Lisbon, Portugal. Zenodo.
##  https://doi.org/10.5281/zenodo.3553442
##
##
## From https://zenodo.org/record/3553442 :
##
##  Vertex attributes:
## 
##     id: unique identifier
##     type: "Person" or "Company"
##     age: years of age for vertices representing a person, "NA" if unknown for a person or if vertex represents a company
##     gender: "Male", "Female" or "NA" (unknown) if vertex represents a person, "NA" if vertex represents a company
##     country: name of country or "NA" (unknown) if vertex represents a company, "NA" if vertex represents a person
##     sector: segment of the economy in which a company operates. "NA" if vertex represents a person or if company's sector is unknown
##     industry: specific business (i.e. subset of sector) in which a company operates. "NA" if vertex represents a person or if company's industry is unknown.
##     employeesnum: number of company's employees. "NA" if vertex represents a person or if company's number of employees is unknown.
## 
##
## Citation for igraph is:
##   
##  Csardi G, Nepusz T (2006). .The igraph software package for complex
##   network research.. InterJournal, Complex Systems, 1695. https://igraph.org.

library(igraph)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  cat("Usage: convertEvtushenkoBoardDirectorAffiliationsToEstimNetDirectedFormat.R board_director.gml.gz\n")
  quit(save="no")
}

##
## network
##
infile <- args[1]

cat("reading ", infile, "...\n")
system.time(g <- read.graph(gzfile(infile), format='gml'))

summary(g)

uniqueIds <- unique(V(g)$id)
numIds <- length(uniqueIds)
stopifnot(min(uniqueIds) == 0)
stopifnot(max(uniqueIds) == numIds-1)
stopifnot(length(V(g)$id) == numIds)

## remove multiple and self edges if any
g <- simplify(g , remove.multiple = TRUE, remove.loops = TRUE)
summary(g)


##
## make sure all Person nodes are first then all Company nodes
##
stopifnot(all(V(g)$type %in% c('Person', 'Company')))
num_Persons <- length(which(V(g)$type == 'Person'))
num_Companies <- length(which(V(g)$type == 'Company'))
cat('num_Persons = ', num_Persons, '\n')
cat('num_Companies = ', num_Companies, '\n')
stopifnot(num_Persons + num_Companies == vcount(g))
stopifnot(all(V(g)$type[1:num_Persons] == 'Person'))

##
## get binary attributes
##

# 0 for Person, 1 for company (duplicates bipartite node type)
# also make a binary gender coded as 0 Male, 1 Female
binattr <- data.frame(company = ifelse(V(g)$type == 'Person', 0, 1),
                      female  = ifelse(V(g)$gender == 'Male', 0,
                                      ifelse(V(g)$gender == 'Female', 1, NA)))
summary(binattr)

##
## get categorical attributes
##
catattr <- data.frame(gender  = ifelse(V(g)$gender == "NA", NA, V(g)$gender),
                      country = ifelse(V(g)$country == "NA", NA, V(g)$country),
                      sector  = ifelse(V(g)$sector == "NA", NA, V(g)$sector),
                      industry= ifelse(V(g)$industry == "NA", NA, V(g)$industry))
print("gender")
catattr$gender <- factor(catattr$gender)
print(levels(catattr$gender))
print("country")
catattr$country <- factor(catattr$country)
print(levels(catattr$country))
print("sector")
catattr$sector <- factor(catattr$sector)
print(levels(catattr$sector))
print("industry")
catattr$industry <- factor(catattr$industry)
print(levels(catattr$industry))

summary(catattr)
catattr$gender <- as.numeric(catattr$gender)
catattr$country <- as.numeric(catattr$country)
catattr$sector <- as.numeric(catattr$sector)
catattr$industry <- as.numeric(catattr$industry)
summary(catattr)
                       

##
## get continuous attributes
##

contattr <- data.frame(age = ifelse(V(g)$age == "NA", NA, as.numeric(V(g)$age)),
                     employeesnum = ifelse(V(g)$employeesnum == "NA", NA, as.numeric(V(g)$employeesnum)))
summary(contattr)

###
### convert type to logical so Pajek bipartite format output works
###
# have to do it this cumbersome way as if do not delete type attr it always
# remains character not logical
V(g)$ltype <- (V(g)$type == 'Company') # TRUE for company, FALSE for person
g <- remove.vertex.attribute(g, 'type')
V(g)$type <- V(g)$ltype
summary(g)
stopifnot(typeof(V(g)$type) == 'logical')
stopifnot(sum(V(g)$type) == num_Companies)

###
### write graph
###


outfilename <- 'evtushenko_directors_bipartite.net'
write.graph(g, outfilename, format="pajek")


##
## write binary attributes
##

write.table(binattr, file = "evtushenko_directors_binattr.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

##
## write categorical attributes
##

write.table(catattr, file = "evtushenko_directors_catattr.txt",
            row.names = FALSE, col.names = TRUE, quote=FALSE)

##
## write continuous attributes
##

write.table(contattr, file = "evtushenko_directors_contattr.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

## end
