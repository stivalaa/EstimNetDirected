#!/usr/bin/env Rscript
##
## File:    convertNorwegianDirectorInterlockDataToEstimNetDirectedFormat.R
## Auithor: Alex Stivala
## Created: April 2024
##
## Read the Seierstad & Opsahl (2011) Norwegian director interlock
## data and convert to format for EstimNetDirected.
##
## Citation for data:
##
## Seierstad, C., & Opsahl, T. (2011). For the few not the many? The
## effects of affirmative action on presence, prominence, and social
## capital of women directors in Norway. Scandinavian journal of
## management, 27(1), 44-54.
##
## Note uses only the network for August 1, 2009, the one used in:
##
## Opsahl, T. (2013). Triadic closure in two-mode networks: Redefining
## the global and local clustering coefficients. Social networks,
## 35(2), 159-167.
##
##
## Output files in cwd (WARNING: overwrites):
##     norwegian_director_interlock_20090801_bipartite.net
##     norwegian_director_interlock_20090801_binattr.txt
##     norwegian_director_interlock_20090801_catattr.txt
##     norwegian_director_interlock_20090801_catattr_strings.txt
##
## Usage:
##    Rscript convertNorwegianDirectorInterlockDataToEstimNetDirectedFormat.R datadir
##    datadir is directory containing the downloaded data from
##    http://www.boardsandgender.com/data.php
##    Note that this website is no longer available, but I
##    downloaded all the data on 29 June 2020. This directory contains
##    the data files:
##        net2m_2009-08-01.txt
##        data_companies.txt
##        data_people.txt
##

library(igraph)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  cat("Usage: convertNorwegianDirectorInterlockDataToEstimNetDirectedFormat.R datadir\n")
  quit(save="no")
}
datadir <- args[1]


##
## read bipatrtite network in tnet format (two column edgelist; note that
## the identifiers for the two modes are not distinct)
##
dat <-read.table(file.path(datadir, 'net2m_2009-08-01.txt'),
                 stringsAsFactors=FALSE, header=FALSE)

## Make string forms of identifiers so (a) the people and company identifiers
## are distinct, and (b) they can be used with igraph::make_bipartite_graph()
## aka igraph::graph.bipartite()

dat$PersonID <- paste("Person", dat$V2, sep='.')
dat$CompanyID <- paste("Company", dat$V1, sep='.')

##
## Make logical node_types vector for make_bipartite_graph().
## It is a named vector where the names are the node names used
## in construting the bipartite graph.
## We will put all the persons (FALSE) first, followed by all the
## companies (TRUE).
## Note: needs igraph version higher than 1.2.11 to allow named nodes and
## types in make_bipartite_graph() (aka graph.bipartite()).
## Written using version 1.3.5.
## Bipartite network node attribute type is FALSE for persons (directors)
## and TRUE for companies (boards).
##
num_Persons <- length(unique(dat$PersonID))
num_Companies <- length(unique(dat$CompanyID))
node_types <- c(rep(FALSE, num_Persons), rep(TRUE, num_Companies))
names(node_types) <- c(unique(dat$PersonID), unique(dat$CompanyID))

##
## Build director affiliation network
##
## Now make the edges vector in the correct format for make_bipartite_graph()
## I.e.:
##  "A vector defining the edges, the first edge points from the first
##   element to the second, the second edge from the third to the fourth, etc.
##   For a numeric vector, these are interpreted as internal vertex ids.
##   For character vectors, they are interpreted as vertex names. "
## [https://search.r-project.org/CRAN/refmans/igraph/html/make_graph.html]
##
## For the transpose and concatenate method to convert the two column
## edge list to this format, see
## https://stackoverflow.com/questions/41051823/convert-a-two-column-matrix-into-a-comma-separated-vector
##
gedges <- c(t(as.matrix(dat[,c("PersonID", "CompanyID")])))
cat('num_Persons = ', num_Persons, '\n')
cat('num_Companies = ', num_Companies, '\n')
## make the bipartite graph
g <- graph.bipartite(types = node_types, edges = gedges, directed = FALSE)
print(summary(g))
stopifnot(typeof(V(g)$type) == 'logical')
stopifnot(sum(V(g)$type) == num_Companies)
stopifnot(num_Persons + num_Companies == vcount(g))
stopifnot(all(V(g)$type[1:num_Persons] == FALSE))
stopifnot(all(V(g)$type[(1+num_Persons):vcount(g)] == TRUE))


# there can be no self edges (since it is bipartite)
stopifnot(!any_loop(g))

## remove multiple edges (these can occur in this data as each row in the
## table is an appointment, so can be appointed to same board in difference
## capacities for example)
print('removing multiple edges...')
g <- simplify(g , remove.multiple = TRUE, remove.loops = FALSE)
summary(g)

##
## Read data on people (directors)
##

## data_people.dat looks like this:
##   $ head /cygdrive/C/Users/alexd/switchdrive/Institution/USI/datasets/interlocking_directorates/Seierstad/data_people.txt
##   "id" "name" "gender"
##   1 "Aage Jakobsen" 1
##   2 "Aage Johan Rem▒y" 1
##   3 "Aage Rasmus Bjelland Figenschou" 1
##   4 "Aagot Irene Skjeldal" 2
##   5 "Aase Gundersen" 2
##   6 "Aase ▒verland" 2
##   7 "Aasmund Fr▒seth" 1
##   8 "Aasmund Rygnestad" 1
##   9 "Aasulv Tveitereid" 1
##
## gender is coded as 1 for male and 2 for female.

data_people <- read.csv(file.path(datadir, 'data_people.txt'),
                        sep=' ', header=TRUE, stringsAsFactors=FALSE)

## Convert identifier to string form to match network data
data_people$PersonID <- paste("Person", data_people$id, sep='.')


##
## Read data on companies
##
## data_companies.txt looks like this:

##   $ head /cygdrive/C/Users/alexd/switchdrive/Institution/USI/datasets/interlocking_directorates/Seierstad/data_companies.txt
##   1       "879447992"     "24SEVENOFFICE ASA"     "0667 OSLO"
##   2       "990031479"     "A-COM NORGE ASA"       "0355 OSLO"
##   3       "890687792"     "ABERDEEN EIENDOMSFOND ASIA ASA"        "0230 OSLO"
##   4       "989761390"     "ABERDEEN EIENDOMSFOND NORDEN/BALTIKUM ASA"     "0255 OSLO"
##   5       "988671258"     "ABERDEEN EIENDOMSFOND NORGE II ASA"    "0255 OSLO"
##   6       "989180797"     "ABERDEEN PROPERTY INVESTORS CORPORATE ASA"     "0230 OSLO"
##   7       "961095026"     "ABG SUNDAL COLLIER HOLDING ASA"        "0250 OSLO"
##   8       "883603362"     "ABG SUNDAL COLLIER NORGE ASA"  "0250 OSLO"
##   9       "989761846"     "ABILITY DRILLING ASA"  "5353 STRAUME"
##   10      "979626878"     "ACTA ASSET MANAGEMENT ASA"     "4006 STAVANGER"

data_companies <- read.csv(file.path(datadir, "data_companies.txt"),
                           sep='\t', header=FALSE, stringsAsFactors=FALSE)
names(data_companies) <- c("id", "organization_number", "fullname", "postcode_and_city")
## Convert identifier to string form to match network data
data_companies$CompanyID <- paste("Company", data_companies$id, sep='.')
## split postcode and city
## everything is always so complicated (and inefficient) in R...
## https://stackoverflow.com/questions/8996134/extract-vectors-from-strsplit-list-without-using-a-loop
## https://stackoverflow.com/questions/8299978/splitting-a-string-on-the-first-space
## replacing first space with colon and split on that as some city names have spaces (but no colons)
data_companies$postcode <- sapply(data_companies$postcode_and_city, function(s) strsplit(sub(' ', ':', s), ':')[[1]], USE.NAMES=FALSE)[1,]
#print(data_companies$postcode)#XXX
data_companies$city <-  sapply(data_companies$postcode_and_city, function(s) strsplit(sub(' ', ':', s), ':')[[1]], USE.NAMES=FALSE)[2,]

#print(data_companies$city)#XXX
#print(data_companies)#XXX

##
## add node attributes
##

person_attrs <- c("gender")
company_attrs <- c("postcode", "city")
all_attrs <- c(person_attrs, company_attrs)

print('adding node attributes to graph...')
# first set all to NA
for (colname in all_attrs) {
  g <- set.vertex.attribute(g, colname, value = NA)
}

##
## check for inconsistent data: make sure there is only one
## entry for each identifier.
##
print('verifying that person attributes are consistent...')
for (attrname in person_attrs) {
  for (person in unique(data_people$PersonID)) {
    attrs <- data_people[which(data_people$PersonID == person), attrname]
    stopifnot(length(attrs) == 1)
  }
}

print('verifying that company attributes are consistent...')
for (attrname in company_attrs) {
  for (company in unique(data_companies$CompanyID)) {
    attrs <- data_companies[which(data_companies$CompanyID == company), attrname]
    stopifnot(length(attrs) == 1)
  }
}


## person attributes
for (colname in person_attrs) {
  for (personid in unique(V(g)[type == FALSE]$name)) {
    ## note the [1] after which() here, to get the first matching row for
    ## that person: it may match several rows, but we checked above that
    ## these all have the same attribute values, so getting the first is safe.
    val <- data_people[which(data_people$PersonID == personid)[1], colname]
    g <- set.vertex.attribute(g, colname, V(g)[personid], val)
  }
}

## company attributes
for (colname in company_attrs) {
##  print(colname)#XXX
  for (companyid in unique(V(g)[type == TRUE]$name)) {
    ## note the [1] after which() here, to get the first matching row for
    ## that company: it may match several rows, but we checked above that
    ## these all have the same attribute values, so getting the first is safe.
#    print(companyid)#XXX
    val <- data_companies[which(data_companies$CompanyID == companyid)[1], colname]
#    print(val)#XXX
    g <- set.vertex.attribute(g, colname, V(g)[companyid], val)
  }
}


##
## get categorical attributes
## Person attributes only: gender
## Company attributes only: postcode, city
##
catattr <- data.frame(gender  = ifelse(V(g)$type == FALSE, V(g)$gender, NA),
                      postcode = ifelse(V(g)$type == TRUE, V(g)$postcode, NA),
                      city = ifelse(V(g)$type == TRUE, V(g)$city, NA)
                     )


summary(catattr)

## print the factor levels to stdout for future reference (codebook)
print("gender")
catattr$gender <- factor(catattr$gender)
print(levels(catattr$gender))
print("postcode")
catattr$postcode <- factor(catattr$postcode)
print(levels(catattr$postcode))
print("city")
catattr$city <- factor(catattr$city)
print(levels(catattr$city))


##
## convert categorical attributes to integer values (as written to stdout above)
##
summary(catattr)
## Also make copy of original categorical attributes data frame
## with strings before conversion to factors, for future referene and
## easier data summary table creation for readability (strings not coded)
catattr_strings <- catattr
catattr_strings$city <- gsub(' ', '_', catattr_strings$city) # replace space with underscore for string output
catattr$gender <- as.numeric(catattr$gender)
catattr$postcode <- as.numeric(catattr$postcode)
catattr$city <- as.numeric(catattr$city)

##
## Make binary attributes
##

binattr <- data.frame(female = ifelse(V(g)$gender == 2, 1, 0))
summary(binattr)

## 
## Also add categorical attributes for male (1 for male, NA otherwise)
## and female (1 for female, NA otherwise).
## so that they can be used as categorical attributes with
## e.g. BipartiteNodematchBetaA in order to test differential homophily
## by having for each gender just one value and for the other NA so it
## never counts for matching (while the gender in question always matches
## itself only).
##

catattr$female <- ifelse(binattr$female == 1, 1, NA)
catattr$male <- ifelse(binattr$female == 0, 1, NA)


summary(catattr)



##
## Writing graph in pajek format uses (as of when this 
## script was written, with igraph version 1.3.5) the internal
## integer node identifiers from 1..N which is just as 
## required for ALAAMEE or EstimNetDirected Pajek network format.
##
write.graph(g, 'norwegian_director_interlock_20090801_bipartite.net',
            format='pajek')

## Write categorical attributes
write.table(catattr, 'norwegian_director_interlock_20090801_catattr.txt',
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(catattr_strings, 'norwegian_director_interlock_20090801_catattr_strings.txt',
            row.names=FALSE, col.names=TRUE, quote=FALSE)


## Write binary attributes
write.table(binattr, 'norwegian_director_interlock_20090801_binattr.txt',
            row.names=FALSE, col.names=TRUE, quote=FALSE)

warnings()
