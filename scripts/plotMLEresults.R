#!/usr/bin/Rscript
#
# File:    plotMLEresults.R
# Author:  Alex Stivala
# Created: September 2013
#
#
#
# Read MLE results (from buildresults.sh) and plot graphs
# of  estimates and standard errors
#
# Usage: Rscript plotMLEresults.R [-a|-s] estimation_results_file.txt
#
# If -a is specified then unconverged samples are not removed, and 
# the filename contains _all after basename
#
# If -s is specified then the standard deviation of the estimated theta
# value is used for the standard error estimate rather than the one
# in the results file computed from the covariance matrix and
# the filename contains _sdtheta after basename.
#
# Output files are PostScript files named
#
#  basename-meanse.eps and
#  basename-boxplot.eps
#  basename-histogram.eps
#  basename-stderrhistogram.eps
#
# where basenanme is basename of snowbal_results_file.txt
#
# e.g.: Rscript plotMLEresults.R estimnetdirected_estimates_n2000_sim.txt
# 

#library(reshape)
#library(doBy)
library(ggplot2)
library(grid)
library(gridExtra)

#zSigma <- 1.96 # number of standard deviations for 95% confidence interval
#zSigma <- 2.58 # number of standard deviations for 99% confidence interval
zSigma <- 2 # nominal 95% CI


options(digits=4) # for printing rmse etc. values on x axis label

args <- commandArgs(trailingOnly=TRUE)
remove_unconverged <- TRUE
use_sd_theta <- FALSE
if (args[1] == "-a") {
  results_filename <- args[2]
  remove_unconverged <- FALSE
} else if  (args[1] == "-s") {
  results_filename <- args[2]
  use_sd_theta <- TRUE
} else {
  results_filename <- args[1]
}
basefilename <- sub("(.+)[.].+", "\\1", basename(results_filename))
if (!remove_unconverged) {
  basefilename <- paste(basefilename, '_all', sep='')
}
if (use_sd_theta) {
  basefilename <- paste(basefilename, '_sdtheta', sep='')
}

  
effects <- c('Arc', 'Reciprocity', 'AinS', 'AoutS', 'AKT-T')

# output effect names, corresponding to above
effect_names <- c('Arc', 'Reciprocity', 'AinStar', 'AoutStar', 'AT-T')

# known true values of effects above (for drawing horizontal line on plot)
true_parameters <- c(-4.0, 4.25, -1.0, -0.5, 1.5)

 
D <- read.table(results_filename, header=TRUE, stringsAsFactors=TRUE)

if (use_sd_theta) {
  D$StdErr <- D$sdEstimate # Use sd(theta) as estimated standard eror
}
  

if (nrow(D[which(D$Effect == 'Receiver'),]) > 0) {
    # if binary attribute present then add the binary attribute effects true values      
    effects <- c('Arc', 'Reciprocity', 'AinS', 'AoutS', 'AKT-T', 'A2P-TD', 'Receiver', 'Sender', 'Interaction')
    effect_names <- c('Arc', 'Reciprocity', 'AinStar', 'AoutStar', 'AKT-T', 'A2P-TD', 'Receiver', 'Sender', 'Interaction')
    true_parameters <- c(-1.0, 4.25, -2.0, -1.5, 0.6, -0.15, 1.0, 1.5, 2.0)
} else if (nrow(D[which(D$Effect == 'Matching'),]) > 0) {
    # if categorical attribute present then add the categorical attribute effects true values      
    effects <- c('Arc', 'Reciprocity', 'AinS', 'AoutS', 'AKT-T', 'A2P-TD', 'Matching', 'MatchingReciprocity')
    effect_names <- c('Arc', 'Reciprocity', 'AinStar', 'AoutStar', 'AKT-T', 'A2P-TD', 'Matching', 'Matching reciprocity')
    true_parameters <- c(-1.0, 4.25, -2.0, -1.5, 1.0, -0.15, 1.5, 2.0)
}

# sort by sampleId so samples arranged horizontally in consistent order
D$sampleId <- as.factor(D$sampleId)
D <- D[sort.list(D$sampleId),]


plotlist <- list()
boxplotlist <- list()
esthistogramlist <- list()
histogramlist <- list()
for (i in 1:length(effects)) {
    effect <- effects[i]
    effect_name <- effect_names[i]
    De <- D[which(D$Effect == effect),]
    print(effect)#XXX

    if (length(De$Effect) == 0) {
        cat("skipping effect ", effect, " not present in data\n")
        # put an empty plot in list at position where this effect would
        # be, to make position of effect plot in grid same, easier to compare
        emptyplot <- ggplot(data.frame()) + geom_blank() +theme(panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.background=element_blank(),strip.background=element_blank(),axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank()) +xlim(0,1)+ylim(0,1) +annotate("text",label=paste("no", effect, "parameter",sep=' '), x=0.5,y=0.5,colour="grey")
        plotlist <- c(plotlist, list(emptyplot))
        boxplotlist <- c(boxplotlist, list(emptyplot))
        histogramlist <- c(histogramlist, list(emptyplot))
        next
    }

    if (remove_unconverged) {
        # remove unconverged samples 
        if ('t.ratio' %in% names(De)) {
            oldn <- nrow(De)
            De <- De[which(De$t.ratio <= 0.1), ]
            if (nrow(De) != oldn) {
              write(paste('removed',oldn-nrow(De),'unconverged estimates\n'),stderr())
            }
        }
    }

        # remove samples with NaN or infinite Estimate
        De <- De[which(!is.na(De$Estimate)), ]
        De <- De[which(!is.infinite(De$Estimate)),]
        De <- De[which(abs(De$Estimate) < 1e03),] # sometimes not inf but still huge
        # remove samples with zero or extreme values as StdErr 
        ## De <- De[which(De$StdErr != 0),]
        ## De <- De[which(!is.na(De$StdErr)),]
        ## De <- De[which(!is.infinite(De$StdErr)),]
        ## De <- De[which(De$StdErr < 1e10),]                                
    
    if (nrow(De) == 0) {
        cat("skipping effect ", effect, ", no converged estimates\n")
        # put an empty plot in list at position where this effect would
        # be, to make position of effect plot in grid same, easier to compare
        emptyplot <- ggplot(data.frame()) + geom_blank() +theme(panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.background=element_blank(),strip.background=element_blank(),axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank()) +xlim(0,1)+ylim(0,1) +annotate("text",label=paste("no", effect, "parameter",sep=' '), x=0.5,y=0.5,colour="grey")
        plotlist <- c(plotlist, list(emptyplot))
        boxplotlist <- c(boxplotlist, list(emptyplot))
        histogramlist <- c(histogramlist, list(emptyplot))
        next
    }

    if (effect == "R_Attribute1"  && substr(basefilename, 1, 7) == "statnet") {
                                        # nodefactor in statnet appears to be not quite the same as R (Activity) in PNet, adjust by subtracting 0.5 from estimate (log-odds) for the statnet value
        De$Estimate <- De$Estimate - 0.5
    }
    
    rmse <- sqrt(mean((De$Estimate - true_parameters[i])^2))
    bias <- mean(De$Estimate - true_parameters[i])

    # count number of times the CI includes true value
    num_in_ci <- sum((De$Estimate < true_parameters[i] & De$Estimate + zSigma*De$StdErr >= true_parameters[i]) | (De$Estimate >= true_parameters[i] & De$Estimate - zSigma*De$StdErr <= true_parameters[i]))
    N <- length(De$Estimate)
    perc_in_ci <- 100 * num_in_ci / length(De$Estimate)

    if (true_parameters[i] == 0.0) {
        # for purely inference (sign and significance, not actual
        # value of estimate) compute False Postivie rate, as the
        # number of times the estimate CI does not include zero
        false_positive_count <- sum(
            (De$Estimate < 0 & De$Estimate + zSigma*De$StdErr < 0) |
            (De$Estimate > 0 & De$Estimate - zSigma*De$StdErr > 0) )
        false_positive_perc <- 100 * false_positive_count / length(De$Estimate)
    } else {
        # for purely inference (sign and significance, not actual
        # value of estimate) compute False Negative rate, as the
        # number of times the estimate is the right sign, but the
        # CI includes zero; or, is the wrong sign.
        false_negative_count <-
            sum( (sign(De$Estimate) == sign(true_parameters[i]) & 
                  ((De$Estimate < 0 & De$Estimate + zSigma*De$StdErr >= 0) |
                   (De$Estimate >= 0 & De$Estimate - zSigma*De$StdErr <= 0)) ) |
                sign(De$Estimate) != sign(true_parameters[i]) )
        false_negative_perc <- 100 * false_negative_count / length(De$Estimate)
    }
    
    # plot each sample estimate with error bar for CI 
    p <- ggplot(De, aes(x = sampleId, y = Estimate))
    p <- p + geom_point()
    p <- p + geom_errorbar(aes(ymax = Estimate + zSigma*StdErr,
                               ymin = Estimate - zSigma*StdErr))
    p <- p + geom_hline(yintercept = true_parameters[i], colour='red')
    p <- p + theme_bw()
    p <- p + theme(panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.background = element_blank(),
                   strip.background = element_blank(),
                   panel.border = element_rect(color = 'black')
                   )
    p <- p + ylab(effect_names[i])
    p <- p + theme(axis.ticks.x = element_blank())
    p <- p + scale_x_discrete(labels=NULL)
    if (true_parameters[i] == 0.0) {
        p <- p + xlab(bquote(atop(list(N[c] == .(N), bias == .(bias), RMSE == .(rmse)), list("% In CI" == .(perc_in_ci), "FPR%" == .(false_positive_perc)))))
    }  else {
        p <- p + xlab(bquote(atop(list(N[c] == .(N), bias == .(bias), RMSE == .(rmse)), list("% In CI" == .(perc_in_ci), "FNR%" == .(false_negative_perc)))))
    }
    plotlist <- c(plotlist, list(p))

    # plot boxplot of all sample estimates
    p <- ggplot(De, aes(x=Effect,y=Estimate))
    p <- p + geom_boxplot()
    p <- p + geom_hline(yintercept = true_parameters[i])
    p <- p + theme_bw()
    p <- p + theme(panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.background = element_blank(),
                   strip.background = element_blank(),
                   panel.border = element_rect(color = 'black')
                   )
    p <- p + ylab(effect_names[i])
    p <- p + theme(axis.ticks.x = element_blank())
    p <- p + scale_x_discrete(labels=NULL)
    if (true_parameters[i] == 0.0) {
        p <- p + xlab(bquote(atop(list(N[c] == .(N), bias == .(bias), RMSE == .(rmse)), list("% In CI" == .(perc_in_ci), "FPR%" == .(false_positive_perc)))))
    } else {
        p <- p + xlab(bquote(atop(list(N[c] == .(N), bias == .(bias), RMSE == .(rmse)), list("% In CI" == .(perc_in_ci), "FNR%" == .(false_negative_perc)))))
    }
    boxplotlist <- c(boxplotlist, list(p))

    # plot histogram of all sample estimates
    p <- ggplot(De, aes(x=Estimate))
    p <- p + geom_histogram(colour='black', fill='white')
    p <- p + ggtitle(effect_names[i])
    p <- p + geom_vline(xintercept = true_parameters[i], colour='red')
    p <- p + geom_vline(xintercept = mean(De$Estimate), colour='blue', linetype='longdash')
    p <- p + theme_bw()
    p <- p + theme(panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.background = element_blank(),
                   strip.background = element_blank(),
                   panel.border = element_rect(color = 'black')
                   )
    if (true_parameters[i] == 0.0) {
        p <- p + xlab(bquote(atop(list(N[c] == .(N), bias == .(bias), RMSE == .(rmse)), list("% In CI" == .(perc_in_ci), "FPR%" == .(false_positive_perc)))))
    } else {
        p <- p + xlab(bquote(atop(list(N[c] == .(N), bias == .(bias), RMSE == .(rmse)), list("% In CI" == .(perc_in_ci), "FNR%" == .(false_negative_perc)))))
    }
    esthistogramlist <- c(esthistogramlist, list(p))

    # plot histogram of standard error values
    p <- ggplot(De, aes(x=StdErr))
    p <- p + geom_histogram(colour='black',fill='white')
    p <- p + ggtitle(effect_names[i])
    p <- p + geom_vline(xintercept = mean(De$StdErr), colour='blue', linetype='longdash')
    p <- p + geom_vline(xintercept = sd(De$Estimate), colour='red')
    p <- p + theme_bw()
    p <- p + theme(panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.background = element_blank(),
                   strip.background = element_blank(),
                   panel.border = element_rect(color = 'black')
                   )
    p <- p + xlab('standard error')
    histogramlist <- c(histogramlist, list(p))
}

postscript(paste(basefilename, '-meanse.eps', sep=''), onefile=FALSE,
           paper="special", horizontal=FALSE, width=9, height=6)
do.call(grid.arrange, plotlist)
dev.off()

postscript(paste(basefilename, '-boxplot.eps', sep=''), onefile=FALSE,
           paper="special", horizontal=FALSE, width=9, height=6)
do.call(grid.arrange, boxplotlist)

postscript(paste(basefilename, '-histogram.eps', sep=''), onefile=FALSE,
           paper="special", horizontal=FALSE, width=9, height=6)
do.call(grid.arrange, esthistogramlist)

postscript(paste(basefilename, '-stderrhistogram.eps', sep=''), onefile=FALSE,
           paper="special", horizontal=FALSE, width=9, height=6)
do.call(grid.arrange, histogramlist)

