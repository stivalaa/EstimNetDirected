#!/usr/bin/Rscript
#
# File:    computeEstimNetDirectedCovariance.R
# Author:  Alex Stivala
# Created: February 2017
#
#
#
# Compute covairance matrix EstimNetDirected output files
#
# Usage: Rscript computeEstimNetDirectedCovariance.R thetaPrefix dzAprefix
#    thetaPrefix is prefix of filenames for theta values 
#    dzAprefix is prefix of filenames for dzA values 
#

library(doBy)
library(reshape2)

#zSigma <- 1.96 # number of standard deviations for 95% confidence interval     
zSigma <- 2.00 # number of standard deviations for nominal 95% confidence interval     

# First iteration number to use, to skip over initial burn-in
firstiter = 20000 # skip first 20000 iterations. FIXME some way to determine properly

# subsample along chains, i.e. take only every x'th 
# iteration for some x, to avoid autocorrelated samples
#period = 50 # take every 50th sample

# Tables have parameter or statistics for each iteration of each run
idvars <- c('run', 't')


args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript plotEstimNetDirectedResults.R thetaPrefix dzAprefix\n")
  quit(save="no")
}
theta_prefix <- args[1]
dzA_prefix <- args[2]

idvars <- c('run', 't')


# get parameter estimates from theta files
theta <- NULL
keptcount <- 0
totalruns <- 0
removed_runs <- NULL
for (thetafile in Sys.glob(paste(theta_prefix, "_[0-9]*[.]txt", sep=''))) {
  run <- as.integer(sub(paste(theta_prefix, "_([0-9]+)[.]txt", sep=''), "\\1",
                        thetafile))
  totalruns <- totalruns + 1
  thetarun <- read.table(thetafile, header=TRUE)
  thetarun$run <- run
  if (any(is.nan(as.matrix(thetarun)))) {
    cat("Removed run", run, "due to NaN\n", file=stderr())
    removed_runs <- c(removed_runs, run)
#  } else if (any(abs(as.matrix(thetarun)) > 1e10)) {
#   cat("Removed run", run, "due to huge values\n", file=stderr())
  } else {
    keptcount <- keptcount + 1
    theta <- rbind(theta, thetarun)
  }
}
cat("Using", keptcount, "of", totalruns, "runs\n", file=stderr())
stopifnot(totalruns - keptcount == length(removed_runs))
theta$run <- as.factor(theta$run)


theta <- theta[which(theta$t > firstiter),]
#theta <- theta[seq(1, nrow(theta), by=period),]

theta <- melt(theta, id=idvars)

# get statistics from dzA files
dzA <- NULL
for (dzAfile in Sys.glob(paste(dzA_prefix, "_[0-9]*[.]txt", sep=''))) {
  run <- as.integer(sub(paste(dzA_prefix, "_([0-9]+)[.]txt", sep=''), "\\1",
                        dzAfile))
  if (!(run %in% removed_runs))  {
    dzArun <- read.table(dzAfile, header=TRUE)
    dzArun$run <- run
    dzA <- rbind(dzA, dzArun)
  } else {
    cat("removed run", run, "from dzA\n", file=stderr())
  }
}
dzA$run <- as.factor(dzA$run)

dzA <- dzA[which(dzA$t > firstiter),]
#dzA <- dzA[seq(1, nrow(dzA), by=period),]

paramnames <- names(dzA)[which(!(names(dzA) %in% idvars))]

dzA <- melt(dzA, id=idvars)

# convert data frame to matrix cols params rows time/run
amatrix <- acast(dzA,  run + t ~ variable  , value.var='value')

# compute covariance matrix 
acov <- cov(amatrix)
acov_inv = solve(acov) # solve(A) is matrix inverse of A
est_stderrs <- sqrt(diag(acov_inv)) 


for (paramname in paramnames) {
    thetav <- theta[which(theta$variable == paramname),]
    thetasum <- summaryBy(value ~ paramname, data=thetav, FUN=c(mean, sd))

    if (paramname != "AcceptanceRate") {
      dzAv <- dzA[which(dzA$variable == paramname),]
      dzAsum <- summaryBy(value ~ paramname, data=dzAv, FUN=c(mean, sd))
      t_ratio <- dzAsum$value.mean / dzAsum$value.sd
    } else {
      t_ratio <- NA
    }

    est_stderr <- est_stderrs[paramname]

    signif <- ''
    if (!is.na(t_ratio) && abs(t_ratio) <= 0.3 && abs(thetasum$value.mean) > zSigma*est_stderr) {
      signif <- '*'
    }

    cat(paramname, thetasum$value.mean, thetasum$value.sd, est_stderr, t_ratio, signif, '\n')
}

cat("TotalRuns", totalruns, "\n")
cat("ConvergedRuns", keptcount, "\n")

