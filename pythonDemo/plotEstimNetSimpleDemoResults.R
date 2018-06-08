#!/usr/bin/Rscript
#
# File:    plotEstimNetSimpleDemoResults.R
# Author:  Alex Stivala
# Created: August 2017
#
# $Id: plotEstimNetSimpleDemoResults.R 16 2017-10-04 03:26:39Z stivalaa $
#
# Plot parmaeter values vs iteration (Algorithm S and then Algorithm EE)
# from the simple Python implementation of the algorithms.
#
#
# Usage: Rscript plotEstimNetSimpleDemoResults.R theta_values.txt  dzA_values.txt 
#                                        
# Input is read fro dzA_values.txt and theta_values.txt written by
# EstimNetSimpleDemo.py
# Output is PostScript files theta_vales.eps and dzA_values.eps
# where theta_values and dzA_values are the base filenames of the 
# supplied .txt files (WARNING: overwritten)
#

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  cat("Usage: plotESimtNetSimpleDemoResults.R theta_Filename.txt dzA_Filename.txt\n")
  quit(save="no")
}
theta_filename <- args[1]
dzA_filename <- args[2]
theta_basefilename <- sub("(.+)[.]txt", "\\1", basename(theta_filename))
dzA_basefilename <- sub("(.+)[.]txt", "\\1", basename(dzA_filename))
theta_outfile <- paste(theta_basefilename, 'eps', sep='.')
dzA_outfile <- paste(dzA_basefilename, 'eps', sep='.')

theta <- read.table(theta_filename, header=T)
dzA <-read.table(dzA_filename, header=T)

postscript(theta_outfile)
par(mfrow=c(3,4))
for (i in names(theta)[2:length(names(theta))]) {
    plot(theta[,"t"], theta[,i], xlab="t", ylab=i)
}

dev.off()

postscript(dzA_outfile)
par(mfrow=c(3,3))
for (i in names(dzA)[2:length(names(dzA))]) {
    plot(dzA[,"t"], dzA[,i], xlab="t", ylab=i, main="dzA")
}
dev.off()




