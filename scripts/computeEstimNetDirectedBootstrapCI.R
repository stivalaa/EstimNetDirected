#!/usr/bin/Rscript
##
## File:    computeEstimNetDirectedBootstrapCI.R
## Author:  Alex Stivala
## Created: January 2019
##
## Read theta parameter MCMC estimates and simulated statistcs from
## EstimNetDirected output and use them to estimate parameter values
## and standard errors using bootstrap procedure.
##
##
## To combine the point estimates and estimated standard errors from
## each run, we use a simple fixed effects meta-analysis, specifically
## inverse-variance weighting (See Ch. 4 of Hartung et al., 2008), to
## combine the estimates so that the overall estimate has the minimum
## variance.
##
##
## Usage: Rscript computeEstimNetDirectedBootstrapCI.R thetaPrefix dzAprefix
##    thetaPrefix is prefix of filenames for theta values 
##    dzAprefix is prefix of filenames for dzA values 
##
##
## References:
##
## Angelo Canty and Brian Ripley (2017). boot: Bootstrap R (S-Plus)
## Functions. R package version 1.3-20.
##
## Davison, A. C. & Hinkley, D. V. (1997) Bootstrap Methods and Their
## Applications. Cambridge University Press, Cambridge. ISBN
## 0-521-57391-2
##
## DiCiccio, T. J. and Efron B. (1996) Bootstrap confidence intervals
## (with Discussion). Statistical Science, 11, 189-228
##
## Efron, B. (1987) More efficient bootstrap computations. Journal of
## the American Statistical Association, 55, 79-89.
##

options(width=9999)  # do not line wrap

library(boot)


t_ratio_threshold <- 0.3 # abs t-ratio must be <= this value for convergence

## First iteration number to use, to skip over initial burn-in
firstiter = 10000 # skip first 10000 iterations. FIXME some way to determine properly


##
## function to obtain mean, in format for use in boot statistic wrapper function
##
## Parameters:
##    estimates - vector of point estimates
##    stderrors - vector of associated std. error estimates
## Return value:
##    pooled estimate (mean)
##
mean_estimate <- function(estimates, stderrors) {
    NumSamples <- length(estimates)
    stopifnot(length(stderrors) == NumSamples) ## even though stderrors not used
    pooledEstimate <- mean(estimates)
    return(pooledEstimate)
}

##
## function to obtain mean estimate, to be used as parameter to boot()
##
mean_bootwrapper <- function(data, indices) {
    estimates <- data[indices,1]
    stderrors <- data[indices,2]
    return( mean_estimate(estimates, stderrors) )
}


##
## Use boot package to compute confidence interval with adjusted bootstrap
## percentile (BCa) method.
##
## Parameters:
##     values - vector of esimated values
##
## Return value:
##    named list with two elements:
##       point_estimate   - pooled estimate (mean)
##       lower - lower value of c.i.
##       upper - upper value of c.i.
##
bootstrap <- function(values) {
    Replicates <- 2000  #FIXME increase (slow) [should use parallel]
    NumSamples <- length(values)
    bootresults <- boot(data = cbind(values, rep(NA, NumSamples)),
                        statistic = mean_bootwrapper,
                        R = Replicates)
    print(bootresults)#XXX
    bootCIresult <- boot.ci(bootresults, type="bca")  # 95% C.I. by default
    print(bootCIresult)#XXX
    return(list(point_estimate = bootresults$t0,
                lower = bootCIresult$bca[4],
                upper = bootCIresult$bca[5]))
}


##############################################################################
### Main
##############################################################################


## Tables have parameter or statistics for each iteration of each run,
## and also acceptance rate, we will have to remove these for the calculations
nonParamVars <- c('run', 't', 'AcceptanceRate')

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript computeEstimNetDirectedBootstrapCI.R thetaPrefix dzAprefix\n")
  quit(save="no")
}
theta_prefix <- args[1]
dzA_prefix <- args[2]



## get parameter estimates from theta files
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
##  } else if (any(abs(as.matrix(thetarun)) > 1e10)) {
##   cat("Removed run", run, "due to huge values\n", file=stderr())
  } else {
    keptcount <- keptcount + 1
    theta <- rbind(theta, thetarun)
  }
}

cat("Using", keptcount, "of", totalruns, "runs\n", file=stderr())
stopifnot(totalruns - keptcount == length(removed_runs))

paramnames <- names(theta)[which(!(names(theta) %in% nonParamVars))]

theta <- theta[which(theta$t > firstiter),]

## get statistics from dzA files
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

dzA <- dzA[which(dzA$t > firstiter),]

if (keptcount < totalruns) {
    ## use R factors to renumber to 1..N and then subtract 1
    ## (since we add one later for the original case of starting at 0)
    fac <- factor(unique(theta$run))
    theta$run <- as.integer(factor(theta$run, levels = levels(fac))) - 1
    dzA$run <- as.integer(factor(dzA$run, levels = levels(fac))) - 1
}

num_runs <- length(unique(theta$run))
stopifnot(length(unique(dzA$run)) == num_runs)
## runs are numbered 0..N-1
stopifnot(min(theta$run) == 0)
stopifnot(max(theta$run) == num_runs-1)
stopifnot(min(dzA$run) == 0)
stopifnot(max(dzA$run) == num_runs-1)

## matrix of theta point estimates, each row is a run (each col a parameter)
theta_estimates <- matrix(nrow = num_runs, ncol = length(paramnames))
colnames(theta_estimates) <- paramnames
## also do this for lower and upper CI bound and t-ratio estimates
ci_lower <- matrix(nrow = num_runs, ncol = length(paramnames))
colnames(ci_lower) <- paramnames
ci_upper <- matrix(nrow = num_runs, ncol = length(paramnames))
colnames(ci_upper) <- paramnames
t_ratios <- matrix(nrow = num_runs, ncol = length(paramnames))
colnames(t_ratios) <- paramnames

for (run in unique(theta$run)) {
    this_theta <- theta[which(theta$run == run), paramnames]
    this_dzA <- dzA[which(dzA$run == run), paramnames]


    theta_sd <- sapply(this_theta, sd)
    

    for (paramname in paramnames) {
        
        bootresult <- bootstrap(this_theta[, paramname])

        stopifnot(bootresult$lower <= bootresult$point_estimate)
        stopifnot(bootresult$upper >= bootresult$point_estimate)
        
        ## runs are numbered from 0 so need to add 1 for R matrix indexing
        theta_estimates[run+1, paramname] <- bootresult$point_estimate
        ci_lower[run+1, paramname] <- bootresult$lower
        ci_upper[run+1, paramname] <- bootresult$upper
        ## estimated t-ratio is mean(dzA)/sd(dzA) for each parameter        
        t_ratios[run+1, paramname] <- mean(this_dzA[,paramname]) /
                                        sd(this_dzA[,paramname])

    }

    ## output estimates for this run
    cat('\nRun ', run, '\n')
    for (paramname in paramnames) {
        signif <- ''
        if ( !is.na(t_ratios[run+1, paramname]) &&
             abs(t_ratios[run+1, paramname]) <= t_ratio_threshold &&
             ( (theta_estimates[run+1, paramname] < 0 &&
                ci_upper[run+1, paramname] < 0) ||
               (theta_estimates[run+1, paramname] >= 0 &&
                ci_lower[run+1, paramname] >= 0) ) ) { 
            signif <- '*'
        }
        cat(paramname, theta_estimates[run+1, paramname], theta_sd[paramname],
            ci_lower[run+1, paramname], ci_upper[run+1, paramname],
            t_ratios[run+1, paramname], signif, '\n')
    }
}





## TODO meta-analysis



cat("TotalRuns", totalruns, "\n")
cat("ConvergedRuns", keptcount, "\n")

