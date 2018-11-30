##!/usr/bin/Rscript
##
## File:    convertSNAPpokecToEstimNetDirectedFormat.R
## Author:  Alex Stivala
## Created: October 2018
##
## Read the network and attributes data for the Pokec social network 
## from SNAP https://snap.stanford.edu/data/soc-Pokec.html
## and convert to Pajek format for EstimNetDirected. See documentation in
## https://snap.stanford.edu/data/soc-pokec-readme.txt
## and source citation:
##
## L. Takac, M. Zabovsky. Data Analysis in Public Social Networks,
## International Scientific Conference & International Workshop
## Present Day Trends of Innovations, May 2012 Lomza, Poland.
## https://snap.stanford.edu/data/soc-pokec.pdf
##
## Reference for SNAP collection of data sets:
## 
##@misc{snapnets,
##  author       = {Jure Leskovec and Andrej Krevl},
##  title        = {{SNAP Datasets}: {Stanford} Large Network Dataset Collection},
##  howpublished = {\url{http://snap.stanford.edu/data}},
##  month        = jun,
##  year         = 2014
##}
##
## Usage:
## 
## Rscript convertSNAPpokecToEstimNetDirectedFormat.R [maxdegree]
##
## If maxdegree is specified then nodes with in or out degree larger
## than maxdegree are removed i.e. the output is the network and
## attributes for the subgraph induced by nodes with in-degree and
## out-degree <= maxdegree.
##
## Input files (in cwd):
##    soc-pokec-profiles.txt.gz
##    soc-pokec-relationships.txt.gz
##
## Output files in cwd (WARNING overwritten):
##     soc-pokec-relationships.txt
##     soc-pokec-binattr.txt
##     soc-pokec-catattr.txt
##     soc-pokec-contattr.txt
##

library(igraph)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 1) {
  cat("Usage: convertSNAPpokecToEstimNetDirectedFormat.R [maxdegree]\n")
  quit(save="no")
}
maxdegree_specified <- FALSE
if (length(args) == 1) {
    maxdegree_specified <- TRUE
    maxdegree <- as.integer(args[1])
}


##
## network
##
infile <- "soc-pokec-relationships.txt.gz"
basefilename <- sub("(.+)[.].+", "\\1", basename(infile))
outfilename <- basefilename

cat("reading ", infile, "...\n")
system.time(edgelist <- read.table(gzfile(infile)))

uniqueIds <- unique(c(edgelist$V1, edgelist$V2))
numIds <- length(uniqueIds)
cat('number of unique ids is ', numIds, '\n')
cat('min is ', min(uniqueIds), ' max is ', max(uniqueIds), '\n')

uniqueIds <- unique(c(edgelist$V1, edgelist$V2))
numIds <- length(uniqueIds)
stopifnot(min(uniqueIds) == 1)
stopifnot(max(uniqueIds) == numIds)

## The Pokec data has nodes numbered 1..N (N = 1632803)
stopifnot(min(uniqueIds) == 1)
stopifnot(numIds == max(uniqueIds))

g <- graph.edgelist(as.matrix(edgelist), directed=TRUE)


summary(g)
## remove multiple and self edges
g <- simplify(g , remove.multiple = TRUE, remove.loops = TRUE)
summary(g)

## https://snap.stanford.edu/data/soc-pokec-readme.txt
stopifnot(vcount(g) == 1632803)
stopifnot(ecount(g) == 30622564)


##
## read attributes
##

pokec_profiles <- 'soc-pokec-profiles.txt.gz'
cat("Reading ", pokec_profiles, "...\n")
system.time( pokec <- read.delim(gzfile(pokec_profiles), header = FALSE,
                                 stringsAsFactors = FALSE) )

## https://snap.stanford.edu/data/soc-pokec-readme.txt
names(pokec) <- c('user_id', 'public', 'completion_percentage',
                    'gender', 'region', 'last_login', 'registration',
                    'AGE', 'body', 'I_am_working_in_field',
                    'spoken_languages', 'hobbies',
                    'I_most_enjoy_good_food', 'pets', 'body_type',
                    'my_eyesight', 'eye_color', 'hair_color',
                    'hair_type', 'completed_level_of_education',
                    'favourite_color', 'relation_to_smoking',
                    'relation_to_alcohol', 'sign_in_zodiac',
                    'on_pokec_i_am_looking_for', 'love_is_for_me',
                    'relation_to_casual_sex', 'my_partner_should_be',
                    'marital_status', 'children',
                    'relation_to_children', 'I_like_movies',
                    'I_like_watching_movie', 'I_like_music',
                    'I_mostly_like_listening_to_music',
                    'the_idea_of_good_evening',
                    'I_like_specialties_from_kitchen', 'fun',
                    'I_am_going_to_concerts', 'my_active_sports',
                    'my_passive_sports', 'profession', 'I_like_books',
                    'life_style', 'music', 'cars', 'politics',
                    'relationships', 'art_culture',
                    'hobbies_interests', 'science_technologies',
                    'computers_internet', 'education', 'sport',
                    'movies', 'travelling', 'health',
                    'companies_brands', 'more')


## Make sure it really does line up with the node ids 1..N
stopifnot(nrow(pokec) == numIds)
stopifnot(min(pokec$user_id) == 1)
stopifnot(max(pokec$user_id) == numIds)

## It does not in fact, so sort by user_id so that it does
pokec <- pokec[order(pokec$user_id), ]
stopifnot(pokec$user_id == 1:nrow(pokec))


## TODO Note there are lots of free-form text fields that we could
## possibly try to parse things from, but it would require a lot of
## manual work to verify/curate to leaving this for later...
## we will only use fields that seem to be well-defined for now.

##
## get binary attributes
##
binattr <- pokec[, c("gender",  # bool, 1 - man
                     "public"   # bool, 1 - all friendships are public
                     )]
## 163 rows have "null" for gender so this converts to NA (with warning)
binattr$gender <- as.numeric(binattr$gender) 
V(g)$gender <- binattr$gender 
V(g)$public <- binattr$public

##
## get categorical attributes
##
catattr <- pokec[, c("gender", "region")]  # also make categorical gender
## https://snap.stanford.edu/data/soc-pokec-readme.txt
## region:
##   string, mostly regions in Slovakia (example: "zilinsky kraj,
##   kysucke nove mesto" means county Zilina, town Kysucke Nove Mesto,
##   Slovakia), some foreign countries (example: "zahranicie, 
##   zahranicie - nemecko" means foreign country Germany (nemecko)),
##   some Czech regions (example: "ceska republika, cz - ostravsky 
##   kraj" means Czech Republic, county Ostrava (ostravsky kraj))
## We just make this a factor, looking at the output written by print
## below, it looks reasonable, but is is only a categorical variable
## allowing us to tell if two users are in the same region or not.
## TODO we could recode this so that we can have different variables
## for being in a different country, major city, etc.
catattr$region <- ifelse(catattr$region == "null", NA, catattr$region)
catattr$region <- factor(catattr$region)
print(levels(catattr$region))

V(g)$region <- catattr$region

##
## get continuous attributes
##
contattr <- pokec[, c("AGE",          # integer, 0 - age attribute not set
                      "registration", # datetime, time at which the
                                        # user registered at the site
                      "completion_percentage" # integer, percentage
                                              # proportion of filled
                                              # values
                      )]
names(contattr)[which(names(contattr) == "AGE")] <- "age"
contattr$age <- as.numeric(contattr$age)
contattr$age <- ifelse(contattr$age == 0, NA, contattr$age)
## convert date of registration to days since January 1, 1970.
## https://www.stat.berkeley.edu/~s133/dates.html
contattr$registration <- as.numeric(as.Date(contattr$registration))

V(g)$age <- contattr$age
V(g)$registration <- contattr$registration
V(g)$completion_percentage <- contattr$completion_percentage

###
### remove hub nodes if specified and write graph
###

if (maxdegree_specified) {
    cat('Will use subgraph induced by nodes with degree <= ', maxdegree, '\n')
        
    g <- induced.subgraph(g, V(g)[which(degree(g, mode='in') < maxdegree &
                                        degree(g, mode='out') < maxdegree)])
    summary(g)
}

write.graph(g, outfilename, format="pajek")


##
## write binary attributes
##

## rename gender to male for binary attribute
binattr <- data.frame(male = V(g)$gender,
                      public = V(g)$public)
summary(binattr$male)
summary(binattr$public)
write.table(binattr, file = "soc-pokec-binattr.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

##
## write categorical attributes
##

## also do categorical version of binary public
catattr <- data.frame(gender = V(g)$gender,
                      region = V(g)$region,
                      publiccat = V(g)$public)
summary(catattr$region)
catattr$region <- as.numeric(catattr$region)
write.table(catattr, file = "soc-pokec-catattr.txt",
            row.names = FALSE, col.names = TRUE, quote=FALSE)

##
## write continuous attributes
##
contattr <- data.frame(age = V(g)$age,
                       registration = V(g)$registration,
                       completion_percentage = V(g)$completion_percentage)
summary(contattr$age)
summary(contattr$registration)
write.table(contattr, file = "soc-pokec-contattr.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

