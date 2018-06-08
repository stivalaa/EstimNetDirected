#!/usr/bin/Rscript
#
# File:    plotTotalTimeHistogram.R
# Author:  Alex Stivala
# Created: January 2014
#
#
#
# Read estimation time results (from buildtimestab.sh) and plot
# histograms of estimation elapsed time
#
# Usage: Rscript plotTotalTimeHistogram.R estimations_time_file.txt
#
#
# Output files are PostScript files named
#
#  basename-totaltimehistogram.eps
#
# where basenanme is basename of _totalestimation_time_file.txt
#
# e.g.: Rscript plotTotalTimeHistogram.R smnet_estimation_times_n5000.txt
# 


             
args <- commandArgs(trailingOnly=TRUE)
results_filename <- args[1]
basefilename <- sub("(.+)[.].+", "\\1", basename(results_filename))


D <- read.table(results_filename, header=TRUE, stringsAsFactors=TRUE)

postscript(paste(basefilename, '-totaltimehistogram.eps', sep=''), onefile=FALSE,
           paper="special", horizontal=FALSE, width=9, height=6)

hist(D$elapsedTime/3600,   xlab = 'Elapsed time (hours)', #breaks=20,
     main = paste('mean total estimation time = ',format(mean(D$elapsedTime/3600), digits=2), 'hours'))
     #summary(D$estimationTime/3600))
         
