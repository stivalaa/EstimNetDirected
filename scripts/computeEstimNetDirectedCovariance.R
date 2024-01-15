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
## Usage: Rscript computeEstimNetDirectedCovariance.R [-v] hetaPrefix dzAprefix
##    thetaPrefix is prefix of filenames for theta values 
##    dzAprefix is prefix of filenames for dzA values 
##
##    -v: verbose mode, output covariance matrices
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
## Vats, D., Flegal, J. M., & Jones, G. L. (2019). Multivariate output
## analysis for Markov chain Monte Carlo. Biometrika, 106(2), 321-337.
## [Preprint available as arXiv:1512.07713]
##
## Vats, D., Flegal, J. M., & Jones, G. L. (2018). Strong consistency
## of multivariate spectral variance estimators in Markov chain Monte
## Carlo. Bernoulli, 24(3), 1860-1909.
##

options(width=9999)  # do not line wrap

library(optparse)
library(mcmcse)

alpha = 0.05  # for 95% confidence interval

t_ratio_threshold <- 0.3 # abs t-ratio must be <= this value for convergence

## First iteration number to use, to skip over initial burn-in
firstiter = 10000 # skip first 10000 iterations. FIXME some way to determine properly

##
## inverse-variance weighting to combine estimates and standard error estimates
## of the runs (See Ch. 4 of Hartung et al., 2008),
##
## Parameters:
##    estimates - vector of point estimates
##    stderrs   - vector of standard error estimates
##
## Return value:
##    Named list with components
##       estimate - inverse-variance weighted mean
##       se       - correpsonding estimated standard error
##
inverse_variance_wm <- function(estimates, stderrs) {
    theta_hat  <- sum(estimates / stderrs^2) / sum(1 / stderrs^2)
    sigma2_hat <- 1 / sum(1 / stderrs^2)
    return(list(estimate = theta_hat,
                se = sqrt(sigma2_hat)))
}


##############################################################################
###
### Main
###
##############################################################################


## Tables have parameter or statistics for each iteration of each run,
## and also acceptance rate, we will have to remove these for the calculations
nonParamVars <- c('run', 't', 'AcceptanceRate')

args <- commandArgs(trailingOnly=TRUE)
option_list <- list(
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Print extra output including covariance matrices")  

  )
parser <- OptionParser(
     usage = "Usage: Rscript computeEximNetDirectedCovariance.R [options] thetaPrefix dzAprefix\n",
     option_list = option_list)
arguments <- parse_args(parser, positional_arguments = 2)
opt <- arguments$options
args <- arguments$args

theta_prefix <- args[1]
dzA_prefix <- args[2]

verbose = opt$verbose


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
  paramnames <- names(thetarun)[which(!(names(thetarun) %in% nonParamVars))]
  if (max(thetarun$t) < firstiter) {
    cat("Removed run", run, "due to not enough iterations\n", file=stderr())
    removed_runs <- c(removed_runs, run)
  } else if (any(is.nan(as.matrix(thetarun[,paramnames])))) {
    cat("Removed run", run, "due to NaN\n", file=stderr())
    removed_runs <- c(removed_runs, run)
  } else if (any(abs(as.matrix(thetarun[,paramnames])) > 1e10)) {
    cat("Removed run", run, "due to huge values\n", file=stderr())
    removed_runs <- c(removed_runs, run)
    ## Otherwise get this error in mcse.multi(): 
    ##   Error in eigen(sig.mat, only.values = TRUE) :
    ##     infinite or missing values in 'x'
    ##   Calls: mcse.multi -> eigen
  } else {
    keptcount <- keptcount + 1
    theta <- rbind(theta, thetarun)
  }
}

stopifnot(totalruns - keptcount == length(removed_runs))


theta <- theta[which(theta$t > firstiter),]

## get statistics from dzA files
dzA <- NULL
for (dzAfile in Sys.glob(paste(dzA_prefix, "_[0-9]*[.]txt", sep=''))) {
  run <- as.integer(sub(paste(dzA_prefix, "_([0-9]+)[.]txt", sep=''), "\\1",
                        dzAfile))
  if (!(run %in% removed_runs))  {
    dzArun <- read.table(dzAfile, header=TRUE)
    amatrix <- as.matrix(dzArun[which(dzArun$t > firstiter), paramnames])
    acov <- cov(amatrix)
    ## http://r.789695.n4.nabble.com/Catching-errors-from-solve-with-near-singular-matrices-td4652794.html
    if (rcond(acov) < .Machine$double.eps)  {
        cat("Removed run ", run, " due to computationally singular covariance matrix (possibly degenerate model)\n", file=stderr())
        removed_runs <- c(removed_runs, run)
        ## also have to remove this run from theta data frame now
        theta <- theta[which(theta$run != run), ]
        keptcount <- keptcount - 1
    } else {
      dzArun$run <- run
      dzA <- rbind(dzA, dzArun)
    }
  } else {
      cat("removed run", run, "from dzA\n", file=stderr())
  }
}

cat("Using", keptcount, "of", totalruns, "runs\n", file=stderr())
stopifnot(totalruns - keptcount == length(removed_runs))

dzA <- dzA[which(dzA$t > firstiter),]

if (keptcount < totalruns) {
    ## use R factors to renumber to 1..N and then subtract 1
    ## (since we add one later for the original case of starting at 0)
    fac <- factor(unique(theta$run))
    theta$run <- as.integer(factor(theta$run, levels = levels(fac))) - 1
    dzA$run <- as.integer(factor(dzA$run, levels = levels(fac))) - 1
}

if (keptcount > 0 ) {
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
  
      ## covariance matrix for MCMC error
      ## mcerror <- mcse.initseq(x = this_theta)
      ## (Not using the conservative initial sequence method (Dai & Jones, 2017)
      ## as we often don't seem to have enough samples for this, unlike the
      ## other methods.)
      ## Instead use the (default) "batch means"
      ## method (Jones 2006; Vats et al., 2017; Vats et al., 2018)
      ## Note the r ("lugsail") parameter was only introduced in
      ## mcmcse R package version 1.4-1 (January 2020); we set it to r=1
      ## to get the "vanilla" estimator so it should behave as it did
      ## before the new version.
      ## See https://cran.r-project.org/web/packages/mcmcse/vignettes/mcmcse_vignette.pdf (January 29, 2020)
      mcerror <- mcse.multi(x = this_theta, r=1, method="bm", size="sqroot")
      est_theta <- mcerror$est    # point estimate (mean)
      Nmcmc <- nrow(this_theta) # number of MCMC samples
      ## mcse.multi returns asymptotic covariance matrix so need to divide
      ## by Nmcmc and take sqrt to get MCMC standard error estimate
      ## (see mcmcse vignette pp.4,8,9): "Note: cov returns an estimate
      ## of \Sigma and not \Sigma/n." (p. 8)
      mcmc_cov <- mcerror$cov / Nmcmc  # covariance matrix
  
      ## covariance matrix for ERGM MLE error
      mcerror_dz <- mcse.multi(x = this_dzA, r=1, method="bm", size="sqroot")
      stopifnot(mcerror_dz$nsim == Nmcmc)
      acov <- mcerror_dz$cov / Nmcmc
      mle_cov = solve(acov) # solve(A) is matrix inverse of A
      
      total_cov <- mcmc_cov + mle_cov
      if (verbose) {
        print(total_cov)
      }
      est_stderr <- sqrt(diag(total_cov))
      names(est_stderr) <- paramnames
  
      theta_sd <- sapply(this_theta, sd)
          
      ## estimated t-ratio is mean(dzA)/sd(dzA) for each parameter
      est_t_ratio <- sapply(this_dzA, FUN = function(v) mean(v)/sd(v))
  
      ## runs are numbered from 0 so need to add 1 for R matrix indexing
      theta_estimates[run+1, ] <- est_theta
      se_estimates[run+1, ] <- est_stderr
      t_ratios[run+1, ] <- est_t_ratio
  
      ## get z-score for alpha (e.g. approx. 1.96 for alpha=0.05)
      zSigma <- qnorm(alpha/2, lower.tail=FALSE)
      
      ## output estimates for this run
      cat('\nRun ', run, '\n')
      for (paramname in paramnames) {
        signif <- ''
          if (!is.na(est_t_ratio[paramname]) &&
              abs(est_t_ratio[paramname]) <= t_ratio_threshold &&
              abs(est_theta[paramname]) > zSigma*est_stderr[paramname]) {
            signif <- '*'
          }
          cat(paramname, est_theta[paramname], theta_sd[paramname],
              est_stderr[paramname], est_t_ratio[paramname], signif, '\n')
      }
    }
}
  
##
## meta-analysis (pooling runs by inverse-variance weighted mean)
##

## get z-score for alpha (e.g. approx. 1.96 for alpha=0.05)
alphaPooled <- alpha
zSigma <- qnorm(alphaPooled/2, lower.tail=FALSE)

cat('\nPooled\n')
if (keptcount > 0) {
  for (paramname in paramnames) {
      pooled_est <- inverse_variance_wm(theta_estimates[, paramname],
                                        se_estimates[, paramname])
  
      ## estimated t-ratio mean(dzA)/sd(dzA) for each parameter,
      ## combining all runs
      est_t_ratio <-  mean(dzA[,paramname])/sd(dzA[,paramname])
  
      ## output pooled estimate for this parameter
      signif <- ''
      if (!is.na(est_t_ratio) &&
          abs(est_t_ratio) <= t_ratio_threshold &&
          abs(pooled_est$estimate) > zSigma*pooled_est$se) {
          signif <- '*'
      }
      cat(paramname, pooled_est$estimate, sd(theta[,paramname]), 
          pooled_est$se, est_t_ratio, signif, '\n')
  }
} else {
  for (paramname in paramnames) {
      cat(paramname, NA, NA,
          NA, NA, '', '\n')
  }
}


cat("TotalRuns", totalruns, "\n")
cat("ConvergedRuns", keptcount, "\n")
##cat("alphaPooled", alphaPooled, "\n")

