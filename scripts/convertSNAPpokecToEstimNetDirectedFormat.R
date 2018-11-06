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
## Rscript convertSNAPpokecToEstimNetDirectedFormat.R
##
## Input files (in cwd):
##    soc-pokec-profiles.txt.gz
##    soc-pokec-relationships.txt.gz
##
## Output files in cwd (WARNING overwritten):
##     soc-pokec-relationships.txt
##
##

library(igraph)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 0) {
  cat("Usage: convertSNAPpokecToEstimNetDirectedFormat.R\n")
  quit(save="no")
}


##
## network
##
infile <- "soc-pokec-relationships.txt.gz"
basefilename <- sub("(.+)[.].+", "\\1", basename(infile))
outfilename <- basefilename


edgelist <- read.table(gzfile(infile))

uniqueIds <- unique(c(edgelist$V1, edgelist$V2))
numIds <- length(uniqueIds)
cat('number of unique ids is ', numIds, '\n')
cat('min is ', min(uniqueIds), ' max is ', max(uniqueIds), '\n')

uniqueIds <- unique(c(edgelist$V1, edgelist$V2))
numIds <- length(uniqueIds)
stopifnot(min(uniqueIds) == 0)
stopifnot(max(uniqueIds) == numIds - 1)

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

write.graph(g, outfilename, format="pajek")


##
## read attributes
##

## https://snap.stanford.edu/data/soc-pokec-readme.txt
names(profile) <- c('user_id', 'public', 'completion_percentage',
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



##
## write binary attributes
##

##
## write categorical attributes
##

##
## write continuous attributes
##
