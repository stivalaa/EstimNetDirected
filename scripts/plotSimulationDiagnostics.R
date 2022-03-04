#!/usr/bin/Rscript
##
## File:    plotSimulationDiagnostics.R
## Author:  Alex Stivala
## Created: November 2019
##
## Derived from plotPNetSimulationDiagnostics.R  (ADS Feb. 2014);
## also plotPNetGofDiagnostics.R (ADS Feb. 2014)
##
## Similarly to the SPSS script genreated by PNet simulation or GoF, plot
## scatterplot to show autocorrelation in samples and histograms of network
## statisics, for use on UNIX version instead of the SPSS script.
##
## Loess smoothed line on scatterplot and mean on histogram are plotted
## in blue dashed lines. On histogram 95% CI is ploted as dotted blue lines.
##
##
## Usage: Rscript plotSimulationDiagnostics.R simulation_stats_output.txt
##                                            [obs_stats.txt]
##
## e.g.: Rscript ../scripts/plotSimulationDiagnostics.R stats_sim_n1000_binattr_sample.txt  obs_stats_n1000_sample_0.txt
##
## If the optional observed stats filename is specified, then the
## observed stats of a single network are read from this and plotted
## in red on the plots for comparison to the simulated stats.
## These observed stats are from the output of EstimNetDirected with
## computeStats = TRUE in file specified by observedStatsFilePrefix.
## The t-ratios of the statistics are also computed and written to stdout.
##
## Output is postscrpt file basename.eps where basename is from the input
## file e.g. stats_sim_n2000_sample-plots.eps
##
##
library(grid)
library(gridExtra)
library(ggplot2)
library(reshape)

zSigma <- 1.96 # number of standard deviations for 95% confidence interval

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1 || length(args) > 2) {
    cat("Usage: Rscript plotSimulationDiagnostics.R simulation_stats.txt [obs_stats.txt]\n")
    quit(save='no')
}
simstats_filename <- args[1]
basefilename <- sub("(.+)[.].+", "\\1", basename(simstats_filename))
do_obs <- FALSE
if (length(args) == 2) {
  obs_stats_filename <- args[2]
  do_obs <- TRUE
}


simstats <- read.table(simstats_filename, header=TRUE, stringsAsFactors=FALSE)

statnames <- names(simstats)[names(simstats) != "t"]

if (do_obs) {
  obsstats <- read.table(obs_stats_filename, header=TRUE, stringsAsFactors=FALSE)
  stopifnot(nrow(obsstats) == 1)
}

simstats <- melt(simstats, id=c('t'))
plotlist <- list()
if (do_obs) {
  cat('effect','observed', 'mean', 'sd', 't_ratio', '\n')
}
for (statname in statnames) {
    simstats_statname <- simstats[which(simstats$variable == statname),]

    p <- ggplot(simstats_statname, aes(x=t, y=value))
    p <- p + geom_point()
    p <- p + geom_smooth(method = loess, color = "blue", linetype = "dashed",
                         se = FALSE)
    if (do_obs && statname != "AcceptanceRate") {
      p <- p + geom_hline(yintercept = obsstats[1,statname],
                          color = "red")
    }
    p <- p + xlab('t')
    p <- p + ylab(statname)
    p <- p + scale_x_continuous(guide = guide_axis(check.overlap = TRUE,
                                                   #n.dodge = 2,
                                                   #angle = 90
                                                  ))
    plotlist <- c(plotlist, list(p))

    p <- ggplot(simstats_statname, aes(x=value))
    p <- p + geom_histogram()
    p <- p + geom_vline(aes(xintercept = mean(value)), color = "blue",
                        linetype = "dashed")
    p <- p + geom_vline(aes(xintercept = mean(value) -
                            zSigma*sd(value)), 
                        colour='blue', linetype='dotted')
    p <- p + geom_vline(aes(xintercept = mean(value) +
                            zSigma*sd(value)), 
                        colour='blue', linetype='dotted')
    if (do_obs && statname != "AcceptanceRate") {
      p <- p + geom_vline(xintercept = obsstats[1,statname],
                          color = "red")
    }
    p <- p + xlab(statname)
    p <- p + scale_x_continuous(guide = guide_axis(check.overlap = TRUE,
                                                   #n.dodge = 2,
                                                   #angle = 90
                                                  ))
    plotlist <- c(plotlist, list(p))

    ## compute t-ratio
    if (do_obs && statname != "AcceptanceRate") {
      simstatvalues<- simstats[which(simstats$variable == statname),"value"]
      t_ratio <- (mean(simstatvalues) - obsstats[1,statname])/sd(simstatvalues)
      cat(statname, obsstats[1,statname], mean(simstatvalues),
           sd(simstatvalues), t_ratio, '\n')
    }
}

postscript(paste(basefilename, '-plots.eps', sep=''), onefile=FALSE,
           paper="special", horizontal=FALSE, width=9, height=6)
do.call(grid.arrange, plotlist)
dev.off()

