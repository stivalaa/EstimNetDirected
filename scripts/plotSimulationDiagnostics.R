#!/usr/bin/Rscript
#
# File:    plotSimulationDiagnostics.R
# Author:  Alex Stivala
# Created: February 2014
#
# Derived from plotPNetSimulationDiagnostics.R  (ADS Feb. 2014).
#
# Similarly to the SPSS script genreated by PNet simulation or GoF, plot
# scatterplot to show autocorrelation in samples and histograms of network
# statisics, for use on UNIX version instead of the SPSS script.
#
#
# Usage: Rscript plotSimulationDiagnostics.R simulation_stats_output.txt
#
# e.g.: Rscript plotSimulationDiagnostics.R stats_sim_n2000_sample.txt
#
# Output is postscrpt file basename.eps where basename is from the input
# file e.g. stats_sim_n2000_sample-plots.eps
#
#
library(grid)
library(gridExtra)
library(ggplot2)
library(reshape)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
    cat("Usage: Rscript plotSimulationDiagnostics.R simulation_stats.txt\n")
    quit(save='no')
}
simstats_filename <- args[1]
basefilename <- sub("(.+)[.].+", "\\1", basename(simstats_filename))


simstats <- read.table(simstats_filename, header=TRUE, stringsAsFactors=FALSE)

statnames <- names(simstats)[names(simstats) != "t"]

simstats <- melt(simstats, id=c('t'))
plotlist <- list()
for (statname in statnames) {
    simstats_statname <- simstats[which(simstats$variable == statname),]

    p <- ggplot(simstats_statname, aes(x=t, y=value))
    p <- p + geom_point()
    p <- p + xlab('t')
    p <- p + ylab(statname)
    plotlist <- c(plotlist, list(p))

    p <- ggplot(simstats_statname, aes(x=value))
    p <- p + geom_histogram()
    p <- p + xlab(statname)
    plotlist <- c(plotlist, list(p))
}

postscript(paste(basefilename, '-plots.eps', sep=''), onefile=FALSE,
           paper="special", horizontal=FALSE, width=9, height=6)
do.call(grid.arrange, plotlist)
dev.off()

