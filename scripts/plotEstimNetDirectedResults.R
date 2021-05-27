#!/usr/bin/Rscript
#
# File:    plotEstimNetDirectedResults.R
# Author:  Alex Stivala
# Created: October 2017
#
#
#
# Plot theta and dzA values vs iteration (Algorithm S and Algorithm EE)
# from EstimNetDirected_mpi output files
#
#
# Usage: Rscript plotEstimNetDirectedResults.R [-p] thetaPrefix dzAprefix
#    thetaPrefix is prefix of filenames for theta values 
#    dzAprefix is prefix of filenames for dzA values 
#      these files have _x.txt appended by EstimNetDirected, where
#      x is task number
#     -p : use PNG not PDF (for very large plots PDF can be huge or even
#          not work correctly with "did not complete the page properly" error
#          when trying to open or convert with ghostscript or PDF viewers)
#          Note uses bitmap() not png() as png() does not work on the Linux
#          cluster I use. See:
#          https://stackoverflow.com/questions/24999983/r-unable-to-start-device-png-capabilities-has-true-for-png
#
# Output files are thetaPrefix.pdf and dzAprefix.pdf
# (or .png if -p option is specified)
# WARNING: output files are overwritten
#
# Example:
#    Rscript plotEstimNetDirectedResults.R theta_sim_sample dzA_sim_sample
#      will use files theta_sim_sample_0.txt dzA_sim_sample_1.txt etc.
#
#
library(grid)
library(gridExtra)
library(ggplot2)
library(reshape)
library(doBy)
library(scales)

#zSigma <- 2.00 # number of standard deviations for 95% confidence interval
zSigma <- 1.96 # number of standard deviations for nominal 95% confidence interval

firstiter = 10000 # XXX skip first 10000 iterations

# http://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
orig_scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}
my_scientific_10 <- function(x) {
# also remove + and leading 0 in exponennt
  parse( text=gsub("e", " %*% 10^", gsub("e[+]0", "e", scientific_format()(x))) )
   
}


args <- commandArgs(trailingOnly=TRUE)
basearg <- 0
use_png <- FALSE
if (length(args) > 3) {
  cat("Usage: Rscript plotEstimNetDirectedResults.R [-p] thetaPrefix dzAprefix\n")
  quit(save="no")
} else if (length(args) == 3) {
  if (args[1] == "-p") {
    use_png <- TRUE
  } else {
  cat("Usage: Rscript plotEstimNetDirectedResults.R [-p] thetaPrefix dzAprefix\n")
  quit(save="no")
  }
  basearg <- basearg + 1
}
theta_prefix <- args[basearg+1]
dzA_prefix <- args[basearg+2]

idvars <- c('run', 't')

if (use_png) {
  theta_outfilename <- paste(theta_prefix, "png", sep='.')
  dzA_outfilename <- paste(dzA_prefix, "png", sep='.')
} else {
  theta_outfilename <- paste(theta_prefix, "pdf", sep='.')
  dzA_outfilename <- paste(dzA_prefix, "pdf", sep='.')
}

theta <- NULL
for (thetafile in Sys.glob(paste(theta_prefix, "_[0-9]*[.]txt", sep=''))) {
  print(thetafile) #XXX
  run <- as.integer(sub(paste(theta_prefix, "_([0-9]+)[.]txt", sep=''), "\\1",
                        thetafile))
  print(run)#XXX
  thetarun <- read.table(thetafile, header=TRUE)
  thetarun$run <- run
  theta <- rbind(theta, thetarun)
}
theta$run <- as.factor(theta$run)
paramnames <- names(theta)[which(!(names(theta) %in% idvars))]
theta <- melt(theta, id=idvars)
plotlist <- list()
for (paramname in paramnames) {
    thetav <- theta[which(theta$variable == paramname),]
    p <- ggplot(thetav, aes(x=t, y=value, colour=run))
    thetasum <- summaryBy(value ~ paramname,
                        data=thetav[which(thetav$t > firstiter), ],
                        FUN=c(mean, sd))

    cat(paramname, thetasum$value.mean, thetasum$value.sd, '\n')
    
    p <- p + geom_point(alpha=1/4)
    p <- p + geom_hline(yintercept=thetasum$value.mean)
    p <- p + geom_hline(yintercept=thetasum$value.mean + zSigma*thetasum$value.sd, linetype='longdash')
    p <- p + geom_hline(yintercept=thetasum$value.mean - zSigma*thetasum$value.sd, linetype='longdash')
    ## uses too much space: p <- p + xlab('t')
    p <- p + ylab(paramname)
#    p <- p + guides(colour = guide_legend(override.aes = list(alpha=1)))
    p <- p + theme(legend.position = "none") # legend with run just makes a mess
    p <- p + theme(axis.title = element_text(size = 8),
                   axis.text = element_text(size = 8),
                   axis.title.x = element_blank(),
                   axis.text.x  = element_text(angle = 45))
    p <- p + scale_x_continuous(labels = my_scientific_10)
    plotlist <- c(plotlist, list(p))
}

print(theta_outfilename)#XXX
if (use_png)  {
  bitmap(theta_outfilename, type='pngalpha', width=9, height=6)
} else {
  # use PDF for transparancy (alpha) not supported by postscript
  pdf(theta_outfilename, onefile=FALSE,
             paper="special", width=9, height=6)
}
do.call(grid.arrange, plotlist)
dev.off()

dzA <- NULL
for (dzAfile in Sys.glob(paste(dzA_prefix, "_[0-9]*[.]txt", sep=''))) {
  print(dzAfile) #XXX
  run <- as.integer(sub(paste(dzA_prefix, "_([0-9]+)[.]txt", sep=''), "\\1",
                        dzAfile))
  print(run)#XXX
  dzArun <- read.table(dzAfile, header=TRUE)
  if (nrow(dzArun) == 0) {
     cat("skipping run ", run, " as no data\n")
     next
   }
  dzArun$run <- run
  dzA <- rbind(dzA, dzArun)
}
dzA$run <- as.factor(dzA$run)
paramnames <- names(dzA)[which(!(names(dzA) %in% idvars))]
dzA <- melt(dzA, id=idvars)
plotlist <- list()
for (paramname in paramnames) {
    dzAv <- dzA[which(dzA$variable == paramname),]
    p <- ggplot(dzAv, aes(x=t, y=value, colour=run))
    dzAsum <- summaryBy(value ~ paramname,
                        data=dzAv[which(dzAv$t > firstiter), ],
                        FUN=c(mean, sd))

    p <- p + geom_point(alpha=1/4)
    p <- p + geom_hline(yintercept=dzAsum$value.mean)
    p <- p + geom_hline(yintercept=dzAsum$value.mean + zSigma*dzAsum$value.sd, linetype='longdash')
    p <- p + geom_hline(yintercept=dzAsum$value.mean - zSigma*dzAsum$value.sd, linetype='longdash')
    ## uses too much space: p <- p + xlab('t')
    p <- p + ylab(paramname)
#    p <- p + guides(colour = guide_legend(override.aes = list(alpha=1)))
    p <- p + theme(legend.position = "none") # legend with run just makes a mess
    p <- p + theme(axis.title = element_text(size = 8),
                   axis.text = element_text(size = 8),
                   axis.title.x = element_blank(),
                   axis.text.x  = element_text(angle = 45))
    p <- p + scale_x_continuous(labels = my_scientific_10)
    plotlist <- c(plotlist, list(p))
}

print(dzA_outfilename)#XXX
if (use_png) {
  bitmap(dzA_outfilename, type='pngalpha', width=9, height=6)
} else {
  # use PDF for transparancy (alpha) not supported by postscript
  pdf(dzA_outfilename, onefile=FALSE,
             paper="special", width=9, height=6)
}
do.call(grid.arrange, plotlist)
dev.off()

