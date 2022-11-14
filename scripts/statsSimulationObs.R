#!/usr/bin/Rscript
##
## File:    statsSimulationsObs.R
## Author:  Alex Stivala
## Created: September 2021
##
## Read sufficient statistics of simulated netwoks from the
## output of SimualteERGM, and sufficient statistics from observed
## network (output of EstimNetDirected for example) and compute
## statistic such as Mahalanobis distance.
##
## Usage: Rscript statsSimulationsObs.R simulation_stats_output.txt
##                                            obs_stats.txt
##
## e.g.: Rscript ../scripts/statsSimulationsObs.R stats_sim_n1000_binattr_sample.txt  obs_stats_n1000_sample_0.txt
##
## The observed stats are from the output of EstimNetDirected with
## computeStats = TRUE in file specified by observedStatsFilePrefix.
##
##
##

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    cat("Usage: Rscript statsSimulationsObs.R simulation_stats.txt obs_stats.txt\n")
    quit(save='no')
}
simstats_filename <- args[1]
basefilename <- sub("(.+)[.].+", "\\1", basename(simstats_filename))
obs_stats_filename <- args[2]

simstats <- read.table(simstats_filename, header=TRUE, stringsAsFactors=FALSE)
statnames <- names(simstats)[!(names(simstats) %in% c("t", "AcceptanceRate"))]
obsstats <- read.table(obs_stats_filename, header=TRUE, stringsAsFactors=FALSE)
stopifnot(nrow(obsstats) == 1)

simstats <- simstats[, names(simstats) %in% statnames] # keep only the stats
simstats <- simstats[, order(names(simstats))] # columns in order by stat name
obsstats <- obsstats[, order(names(obsstats))] # columns in order by stat name
#print(ncol(simstats))#XXX
#print(ncol(obsstats))#XXX
#print(statnames)#XXX
#head(simstats)#XXX
#print(obsstats)#XXX
stopifnot(all(names(simstats) == names(obsstats)))

## convert to matrix where columns are the statistics 
## and rows are simulated networks
## withsingle row for the single observation as last row
stats_matrix <- as.matrix(simstats)
#print(stats_matrix)#XXX
stats_matrix <- rbind(stats_matrix, as.matrix(obsstats))
#print(stats_matrix)#XXX

covstats <- cov(stats_matrix)
inverted_cov_stats_matrix <- solve(covstats) # inverse of covstats
## mahalanobis() returns squred Mahalanobis distance so do sqrt()
mdist <- sqrt(mahalanobis(stats_matrix, colMeans(stats_matrix), inverted_cov_stats_matrix, inverted=TRUE))
#print(mdist)#XXX

obs_mdist <- mdist[length(mdist)] # Mahalanobis distance of observed (last in vector from mahalanobis())

cat("statsistics = ", statnames, "\n")
cat("Mahalanobis = ", obs_mdist, "\n")


