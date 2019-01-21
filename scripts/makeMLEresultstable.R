#!/usr/bin/Rscript
#
# File:    makeMLEresultstable.R
# Author:  Alex Stivala
# Created: November 2013
#
#
#
# Read estimation results (from buildresults.sh) and plot graphs
# of estimates and standard errors
#
# Usage: Rscript [-s] makeMLEresultstable.R 
#
#
# Output is to stdout
#
# e.g.: Rscript makeMLEresultstable.R
# 

library(PropCIs) # for Wilson score test for false negative rate

#zSigma <- 1.96 # number of standard deviations for 95% confidence interval
#zSigma <- 2.58 # number of standard deviations for 99% confidence interval
zSigma <- 2 # nominal 95% CI


options(digits=4) # for printing rmse etc. values


results_filenames <- c(#'estimnetdirected_estimates_n2000_sim.txt',
                       'estimnetdirected_estimates_n2000_binattr.txt' 
                       #'estimnetdirected_estimates_n2000_cat3.txt',
                       #'estimnetdirected_estimates_n5000_sim.txt'
                      )

args <- commandArgs(trailingOnly=TRUE)
use_sd_theta <- FALSE
if (length(args) > 0 ) {
   if (args[1] == "-s") {
      use_sd_theta <- TRUE
    }
}  

# now only write the rows, use shell script to sort and add header

for (results_filename in results_filenames) {
  Dorig <- read.table(results_filename, header=TRUE, stringsAsFactors=TRUE)

    # known true values of effects above (for drawing horizontal line on plot)
    true_parameters <- c(-4.0, 4.25, -1.0, -0.5, 1.5)

    effects <- c('Arc', 'Reciprocity', 'AinS', 'AoutS', 'AKT-T')

    effect_names <- c('Arc', 'Reciprocity', 'AinStar', 'AoutStar', 'AT-T')
    num_seedsets <- NA
                
    # get description of data from filename
    error_analysis <- '*UNKNOWN*' # should always be replaced
    fixeddensity <- 'N'
    network_N_list <- NULL # should always be replaced
    attribute_descr <- 'None'
    if (length(grep('fixdensity', results_filename)) > 0) {
      fixeddensity <- 'Y'
    }
    if (length(grep('bootstrapanalysis', results_filename)) > 0) {
      error_analysis <- 'bootstrap'
    }
    else if (length(grep('metaanalysis', results_filename)) > 0 ){
      error_analysis <- 'WLS'
    }
    else if (substr(results_filename, 1, 4) == 'pnet') {
        # not really an error analysis method, but abuse it for PNet esimations
        error_analysis <- 'PNet'
    }
    else if (substr(results_filename, 1, 13) == 'statnet_mcmle') {
        # not really an error analysis method, but abuse it for statnet esimations
        error_analysis <- 'statnet MCMLE'
    }
    else if (substr(results_filename, 1, 16) == 'statnet_stepping') {
        # not really an error analysis method, but abuse it for statnet esimations
        error_analysis <- 'statnet Stepping'
    }
    else if (substr(results_filename, 1, 5) == 'smnet') {
        # not really an error analysis method, but abuse it for SMNet esimations
        error_analysis <- 'SMNet'
    }
    else if (substr(results_filename, 1, 16) == 'estimnetdirected') {
        # not really an error analysis method, but abuse it for EstimNetDirected estimations
        error_analysis <- 'EstimNetDirected'
    }
    else if (any(grepl('one_large_sample', results_filename))) {
        # not really an error analysis method, but abuse it for driect contidional esimations with no pooling (metaanalsysi / bootstrap)
        error_analysis <- 'no pooling'
    }                
    if (length(grep('n5000', results_filename)) > 0 ) {
      network_N_list <- 5000
    } else if (length(grep('n500_', results_filename)) > 0 ) {
      network_N_list <- 500
    } else if (length(grep('n10k', results_filename)) > 0) {
      network_N_list <- 10000
    } else if (length(grep('n1000_', results_filename)) > 0) {
      network_N_list <- 1000
    } else if (any(grepl('n2000_', results_filename, fixed=TRUE))) {
      network_N_list <- 2000
    } else if (any(grepl('n100_', results_filename, fixed=TRUE))) {
      network_N_list <- 100
    } else if (any(grepl('manyN_', results_filename, fixed=TRUE))) {
      network_N_list <- unique(Dorig$nodeCount)
    }

    if (length(grep('binattr', results_filename, fixed=TRUE)) > 0) {
        # if binary attribute present then add the binary attribute effects true values      
        effects <- c('Arc', 'Reciprocity', 'AinS', 'AoutS', 'AKT-T', 'A2P-TD', 'Receiver', 'Sender', 'Interaction')
        effect_names <- c('Arc', 'Reciprocity', 'AinStar', 'AoutStar', 'AT-T', 'A2P-TD', 'Receiver', 'Sender', 'Interaction reciprocity')
        true_parameters <- c(-1.0, 4.25, -2.0, -1.5, 0.6, -0.15, 1.0, 1.5, 2.0)
    } else if (length(grep('cat3', results_filename, fixed=TRUE)) > 0) {
        # if categorical attribute present then add the categorical attribute effects true values      
        effects <- c('Arc', 'Reciprocity', 'AinS', 'AoutS', 'AKT-T', 'A2P-TD', 'Matching', 'MatchingReciprocity')
        effect_names <- c('Arc', 'Reciprocity', 'AinStar', 'AoutStar', 'AKT-T', 'A2P-TD', 'Matching', 'Matching reciprocity')
        true_parameters <- c(-1.0, 4.25, -2.0, -1.5, 1.0, -0.15, 1.5, 2.0)
    }

     
    for (network_N in network_N_list) {
      if ('nodeCount' %in% names(Dorig)) {
          D <- Dorig[which(Dorig$nodeCount == network_N), ]
      }  else {
          D <- Dorig
      }

      if (use_sd_theta) {
         D$StdErr <- D$sdEstimate # Use sd(theta) as estimated standard eror
      }

	    if (length(grep('binattr', results_filename)) > 0) {
	      attribute_descr <- 'Binary' 
	    }
	    else if (length(grep('cat3', results_filename)) > 0) {
        attribute_descr <- 'Categorical'
	    }
	    else if (length(grep('cont', results_filename)) > 0) {
        attribute_descr <- 'Continuous'
	    }


      this_effects <- effects
	    this_true_parameters <- true_parameters

			for (i in 1:length(this_effects)) {
        effect <- this_effects[i]
        De <- D[which(D$Effect == effect),]
        if (length(De$Effect) == 0) {
            next
        }


        # remove unconverged samples 
        if ('t.ratio' %in% names(De)) {
            oldn <- nrow(De)
            De <- De[which(De$t.ratio <= 0.1), ]
            if (nrow(De) != oldn) {
                write(paste('removed',oldn-nrow(De),'unconverged estimates\n'),stderr())
            }
        }

        # remove samples with NaN or infinite Estimate
        De <- De[which(!is.na(De$Estimate)), ]
        De <- De[which(!is.infinite(De$Estimate)),]
        ##De <- De[which(abs(De$Estimate) < 1e03),] # sometimes not inf but still huge
        # remove samples with zero or extreme values as StdErr 
        ## De <- De[which(De$StdErr != 0),]
        ## De <- De[which(!is.na(De$StdErr)),]
        ## De <- De[which(!is.infinite(De$StdErr)),]
        ## De <- De[which(De$StdErr < 1e10),]


        totalRuns <- unique(De$totalRuns)
        stopifnot(length(totalRuns) == 1)
        
        if (effect == "R_Attribute1"  && substr(error_analysis, 1, 7) == "statnet") {
            # nodefactor in statnet appears to be not quite the same as R (Activity) in PNet, adjust by subtracting 0.5 from estimate (log-odds) for the statnet value
            De$Estimate <- De$Estimate - 0.5
        }

        rmse <- sqrt(mean((De$Estimate - this_true_parameters[i])^2))
        bias <- mean(De$Estimate - this_true_parameters[i])

        # count number of times the CI  includes true value
        num_in_ci <- sum((De$Estimate < this_true_parameters[i] & De$Estimate + zSigma*De$StdErr >= this_true_parameters[i]) | (De$Estimate >= this_true_parameters[i] & De$Estimate - zSigma*De$StdErr <= this_true_parameters[i]))
        perc_in_ci <- 100 * num_in_ci / length(De$Estimate)

        # for purely inference (sign and significance, not actual value of estimate)
        # compute False Negative rate, as the number of times the estimate
        # is the right sign, but the CI includes zero; or, is the wrong sign.
            
        false_negative_count <-
            sum( (sign(De$Estimate) == sign(this_true_parameters[i]) & 
                  ((De$Estimate < 0 & De$Estimate + zSigma*De$StdErr >= 0) |
                    (De$Estimate >= 0 & De$Estimate - zSigma*De$StdErr <= 0)) ) |
                  sign(De$Estimate) != sign(this_true_parameters[i]) )
        false_negative_perc <- 100 * false_negative_count / length(De$Estimate)
        confint <- scoreci(false_negative_count, length(De$Estimate), 0.95)
        fnp_lower <- confint$conf.int[1] * 100
        fnp_upper <- confint$conf.int[2] * 100
        cat(network_N, attribute_descr, error_analysis, fixeddensity, 
        gsub('_', ' ', effect), bias, rmse, false_negative_perc,
        mean(De$StdErr), sd(De$Estimate), mean(De$convergedRuns),
                fnp_lower, fnp_upper, num_seedsets,
                length(De$Estimate), totalRuns, perc_in_ci,
        sep=' & ')
        cat('\\\\\n')
        }

    }
}

