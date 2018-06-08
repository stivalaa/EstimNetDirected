#!/usr/bin/Rscript
#
# File:    plotStdErrSdThetaScatterplot.R
# Author:  Alex Stivala
# Created: November 2017
#
#                                        
# Usage: Rscript plotStdErrSdThetaScatterplot.R estimation_results_file.txt
#
#
# Make scatter plot of estimated standard error against
# standard deviation of parameter estimate in each estimation result.
#
# Output file is  PostScript file named
#
#  basename-stderrscatterplot.eps and
#
# where basenanme is basename of estimation_results_file.txt
#
# 
args <- commandArgs(trailingOnly=TRUE)
results_filename <- args[1]
basefilename <- sub("(.+)[.].+", "\\1", basename(results_filename))

D <- read.table(results_filename, header=TRUE)

D$Effect <- as.factor(D$Effect)

postscript(paste(basefilename, '-stderrscatterplot.eps', sep=''),
           onefile=FALSE, paper="special", horizontal=FALSE, width=9, height=9)

plot(D$sdEstimate, D$StdErr,
     col=D$Effect,
     xlab=expression(sd(theta)), ylab='estimated std. error')

legend('topleft', bty='n',
       legend=levels(D$Effect),
       col=seq(nlevels(D$Effect)),
       pch=rep(1, nlevels(D$Effect)))

dev.off()


