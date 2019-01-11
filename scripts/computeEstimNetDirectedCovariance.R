#!/usr/bin/Rscript
##
## File:    computeEstimNetDirectedCovariance.R
## Author:  Alex Stivala
## Created: January 2019
##
## Read theta parameter MCMC estimates and simulated statistcs from
## EstimNetDirected output and use them to estimate parameter values
## and standard errors.
##
## For each independent (parallel) run, the parameter value and MCMC
## covariance are estimated (using the mcmcse package). However as
## well as error due to MCMC, we also have the covariance due to use
## of MLE; this is estimated from the Fisher information matrix as per
## Snijders (2002) [see also Hunter & Handcock (2006)]; the inverse of
## this covariance matrix is then the covariance due to ERGM MLE of
## the parameter esetimate.  These two sources of uncertainty are
## combined by adding the covariance matrices to get the total
## uncertainty in our theta estimates.  From this we then get a
## standard error estimate for that run.
##
## (Note we do not combine all the runs together first, as the MCMC error
## estimation is using the fact that the run is a chain to adjust error
## according to sample size, autocorrelation, etc.).
##
## To combine the point estimates and estimated standard errors from
## each run, we use a simple fixed effects meta-analysis, specifically
## inverse-variance weighting (See Ch. 4 of Hartung et al., 2008), to
## combine the estimates so that the overall estimate has the minimum
## variance.
##
##
## Usage: Rscript computeEstimNetDirectedCovariance.R thetaPrefix dzAprefix
##    thetaPrefix is prefix of filenames for theta values 
##    dzAprefix is prefix of filenames for dzA values 
##
##
## References:
##
## Dai, N., & Jones, G. L. (2017). Multivariate initial sequence
## estimators in Markov chain Monte Carlo. Journal of Multivariate
## Analysis, 159, 184-199.
##
## James M. Flegal, John Hughes, Dootika Vats, and Ning
## Dai. (2017). mcmcse: Monte Carlo Standard Errors for MCMC. R package
## version 1.3-2. Riverside, CA, Denver, CO, Coventry, UK, and
## Minneapolis, MN.
##
## Flegal, J. M., & Jones, G. L. (2010). Batch means and spectral
## variance estimators in Markov chain Monte Carlo. The Annals of
## Statistics, 38(2), 1034-1070.
##
## Hartung, J., Knapp, G., & Sinha, B. K. (2008). Statistical
## meta-analysis with applications. John Wiley & Sons. Hoboken, NJ.
##
## Hunter & Handcock (2006) "Inference in Curved Exponential Family
## Models for Networks" J. Comp. Graph. Stat. 15:3, 565-583
##
## Jones, G. L., Haran, M., Caffo, B. S., & Neath,
## R. (2006). Fixed-width output analysis for Markov chain Monte
## Carlo. Journal of the American Statistical Association, 101(476),
## 1537-1547.
##
## Snijders (2002) "Markov chain Monte Carlo estimation of exponential
## random graph models" J. Social Structure 3(2):1-40).
##
## Vats, D., Flegal, J. M., & Jones, G. L. (2017). Multivariate output
## analysis for Markov chain Monte Carlo. arXiv preprint
## arXiv:1512.07713.
##
## Vats, D., Flegal, J. M., & Jones, G. L. (2018). Strong consistency
## of multivariate spectral variance estimators in Markov chain Monte
## Carlo. Bernoulli, 24(3), 1860-1909.
##

options(width=9999)  # do not line wrap

library(mcmcse)

##zSigma <- 1.96 # number of standard deviations for 95% confidence interval     
zSigma <- 2.00 # number of standard deviations for nominal 95% confidence interval     

t_ratio_threshold <- 0.3 # abs t-ratio must be <= this value for convergence

## First iteration number to use, to skip over initial burn-in
firstiter = 10000 # skip first 10000 iterations. FIXME some way to determine properly


## Tables have parameter or statistics for each iteration of each run,
## and also acceptance rate, we will have to remove these for the calculations
nonParamVars <- c('run', 't', 'AcceptanceRate')

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript computeEstimNetDirectedCovariance.R thetaPrefix dzAprefix\n")
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
## also do this for standard error and t-ratio estimates
se_estimates <- matrix(nrow = num_runs, ncol = length(paramnames))
colnames(se_estimates) <- paramnames
t_ratios <- matrix(nrow = num_runs, ncol = length(paramnames))
colnames(t_ratios) <- paramnames

for (run in unique(theta$run)) {
    this_theta <- theta[which(theta$run == run), paramnames]
    this_dzA <- dzA[which(dzA$run == run), paramnames]

    ## covariance matrix for ERGM MLE error
    acov <- cov(as.matrix(this_dzA))
    mle_cov = solve(acov) # solve(A) is matrix inverse of A

    ## covariance matrix for MCMC error
    ## mcerror <- mcse.initseq(x = this_theta)
    ## (Not using the conservative initial sequence method (Dai & Jones, 2017)
    ## as we often don't seem to have enough samples for this, unlike the
    ## other methods.)
    ## Instead use the (default) "batch means"
    ## method (Jones 2006; Vats et al., 2017; Vats et al., 2018)
    mcerror <- mcse.multi(x = this_theta, method="bm") 
    est_theta <- mcerror$est    # point estimate (mean)
    Nmcmc <- nrow(this_theta) # number of MCMC samples
    ## mcse.multi returns asymptotic covariance matrix so need to divide
    ## by Nmcmc and take sqrt to get MCMC standard error estimate
    ## (see mcmcse vignette pp.4,8,9): "Note: cov returns an estimate
    ## of \Sigma and not \Sigma/n." (p. 8)
    mcmc_cov <- mcerror$cov / Nmcmc  # covariance matrix

    total_cov <- mcmc_cov + mle_cov
    est_stderr <- sqrt(diag(total_cov))

    # estimated t-ratio is mean(dzA)/sd(dzA) for each parameter
    est_t_ratio <- sapply(this_dzA, FUN = function(v) mean(v)/sd(v))

    # runs are numbered from 0 so need to add 1 for R matrix indexing
    theta_estimates[run+1, ] <- est_theta
    se_estimates[run+1, ] <- est_stderr
    t_ratios[run+1, ] <- est_t_ratio


    cat('\nRun ', run, '\n')
    for (paramname in paramnames) {
        signif <- ''
        if (!is.na(est_t_ratio[paramname]) &&
            abs(est_t_ratio[paramname]) <= t_ratio_threshold &&
            abs(est_theta[paramname]) > zSigma*est_stderr[paramname]) {
            signif <- '*'
        }
        cat(paramname, est_theta[paramname], est_stderr[paramname],
            est_t_ratio[paramname], signif, '\n')
    }
}





## TODO meta-analysis



cat("TotalRuns", totalruns, "\n")
cat("ConvergedRuns", keptcount, "\n")

