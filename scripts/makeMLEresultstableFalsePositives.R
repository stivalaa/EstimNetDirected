#!/usr/bin/Rscript
#
# File:    makeMLEresultstableFalsePositives.R
# Author:  Alex Stivala
# Created: December 2013
#
#
#
# Read snowball sampling results (from buildresults.sh) and make table
# of snowball estimates and standard errors.
# This version, unlike makeMLEresultstable.R, handles the simulated
# networks each with one effect set to zero, to identifify false positvies
# in the inference (rather than false negatives).
#
# Usage: Rscript makeMLEresultstableFalsePositives.R 
#
#
# Output is to stdout
#
# e.g.: Rscript makeMLEresultstableFalsePositives.R
# 

library(PropCIs) # for Wilson score test for false negative rate


#zSigma <- 1.96 # number of standard deviations for 95% confidence interval
#zSigma <- 2.58 # number of standard deviations for 99% confidence interval
zSigma <- 2 # nominal 95% CI


options(digits=4) # for printing rmse etc. values


results_filenames <- c('estimnetdirected_estimates_n2000_binattr_A2P0.txt', 
                       'estimnetdirected_estimates_n2000_binattr_AT0.txt', 
                       'estimnetdirected_estimates_n2000_binattr_interaction0.txt', 
                       'estimnetdirected_estimates_n2000_binattr_reciprocity0.txt', 
                       'estimnetdirected_estimates_n2000_binattr_sender0.txt', 
                       'estimnetdirected_estimates_n2000_binattr_receiver0.txt',
                       'estimnetdirected_estimates_n2000_cat3_reciprocity0.txt',
                       'estimnetdirected_estimates_n2000_cat3_matchingreciprocity0.txt')

 
args <- commandArgs(trailingOnly=TRUE)
use_sd_theta <- FALSE
if (length(args) > 0 ) {
   if (args[1] == "-s") {
      use_sd_theta <- TRUE
    }
}  

# now only write the rows, use shell script to sort and add header

for (results_filename in results_filenames) {
     # known true values of effects above (for drawing horizontal line on plot)
    true_parameters <- c(-4.0, 4.25, -1.0, -0.5, 1.5)

    effects <- c('Arc', 'Reciprocity', 'AinS', 'AoutS', 'AKT-T')

    effect_names <- c('Arc', 'Reciprocity', 'AinStar', 'AoutStar', 'AT-T')
    num_seedsets <- NA


    D <- read.table(results_filename, header=TRUE, stringsAsFactors=TRUE)

    if (use_sd_theta) {
      D$StdErr <- D$sdEstimate # Use sd(theta) as estimated standard eror
    }
    # get description of data from filename
    error_analysis <- '*UNKNOWN*' # should always be replaced
    fixeddensity <- 'N'
    network_N <- NULL # should always be replaced
    attribute_descr <- 'None'


    if (length(grep('n5000', results_filename)) > 0 ) {
      network_N <- 5000
    } else if (length(grep('n500_', results_filename)) > 0 ) {
      network_N <- 500
    } else if (length(grep('n10k', results_filename)) > 0) {
      network_N <- 10000
    } else if (length(grep('n1000_', results_filename)) > 0) {
      network_N <- 1000
    } else if (any(grepl('n2000_', results_filename, fixed=TRUE))) {
      network_N <- 2000
    } else if (any(grepl('n100_', results_filename, fixed=TRUE))) {
      network_N <- 100
    }
    
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

    if (length(grep('binattr', results_filename, fixed=TRUE)) > 0) {
        # if binary attribute present then add the binary attribute effects true values      
      attribute_descr <- 'Binary'
      effects <- c('Arc', 'Reciprocity', 'AinS', 'AoutS', 'AKT-T', 'A2P-TD', 'Receiver', 'Sender', 'Interaction')
      effect_names <- c('Arc', 'Reciprocity', 'AinStar', 'AoutStar', 'AT-T', 'A2P-TD', 'Receiver', 'Sender', 'Interaction reciprocity')
      true_parameters <- c(-1.0, 4.25, -2.0, -1.5, 0.6, -0.15, 1.0, 1.5, 2.0)

      zero_effect <- "*UNKNOWN*" # should always be replced with a valid name
                
      # use filename to see if data had one of the parameters set to zero
      if (length(grep("_reciprocity0", results_filename)) > 0) {
        true_parameters[2] <- 0.0
        zero_effect <- effects[2]
      }
        
      if (length(grep("_AT0", results_filename)) > 0) {
        true_parameters[5] <- 0.0
        zero_effect <- effects[5]
      }
      if (length(grep("_A2P0", results_filename)) > 0) {
        true_parameters[6] <- 0.0
        zero_effect <- effects[6]
      }
      if (length(grep("_receiver0", results_filename)) > 0 ) {
        true_parameters[7] <- 0.0
        zero_effect <- effects[7]
      }
      if (length(grep("_sender0", results_filename)) > 0 ) {
        true_parameters[8] <- 0.0
        zero_effect <- effects[8]
      }
      if (length(grep("_interaction0", results_filename)) > 0 ) {
        true_parameters[9] <- 0.0
        zero_effect <- effects[9]
      }
    } else if (length(grep('cat3', results_filename, fixed=TRUE)) > 0) {
      # if categorical attribute present then add the categorical attribute effects true values      
      attribute_descr <- 'Categorical'
      effects <- c('Arc', 'Reciprocity', 'AinS', 'AoutS', 'AKT-T', 'A2P-TD', 'Matching', 'MatchingReciprocity')
      effect_names <- c('Arc', 'Reciprocity', 'AinStar', 'AoutStar', 'AKT-T', 'A2P-TD', 'Matching', 'Matching reciprocity')
      true_parameters <- c(-1.0, 4.25, -2.0, -1.5, 1.0, -0.15, 1.5, 2.0)
      # use filename to see if data had one of the parameters set to zero
      if (length(grep("_reciprocity0", results_filename)) > 0) {
        true_parameters[2] <- 0.0
        zero_effect <- effects[2]
      }
        
      if (length(grep("_AT0", results_filename)) > 0) {
        true_parameters[5] <- 0.0
        zero_effect <- effects[5]
      }
      if (length(grep("_A2P0", results_filename)) > 0) {
        true_parameters[6] <- 0.0
        zero_effect <- effects[6]
      }
      if (length(grep("_matching0", results_filename)) > 0 ) {
        true_parameters[7] <- 0.0
        zero_effect <- effects[7]
      }
      if (length(grep("_matchingreciprocity0", results_filename)) > 0 ) {
        true_parameters[8] <- 0.0
        zero_effect <- effects[8]
      }
    }
    
    for (i in 1:length(effects)) {
                    effect <- effects[i]
                    De <- D[which(D$Effect == effect),]

                    if (length(De$Effect) == 0) {
                                    next
                    }

                    # remove samples with NaN or infinite Estimate
                    De <- De[which(!is.na(De$Estimate)), ]
                    De <- De[which(!is.infinite(De$Estimate)),]
                    De <- De[which(abs(De$Estimate) < 1e03 ),] # sometimes not inf but still huge
                    # remove samples with zero or extreme values as StdErr 
                    ## De <- De[which(De$StdErr != 0),]
                    ## De <- De[which(!is.na(De$StdErr)),]
                    ## De <- De[which(!is.infinite(De$StdErr)),]
                    ## De <- De[which(De$StdErr < 1e10),]

                    if (length(De$Effect) == 0) {
                        cat("Removed all estimates for", effect,"\n",file=stderr())
                        cat(network_N, attribute_descr, error_analysis, fixeddensity, 
                            gsub('_', ' ', effect), sub('_', ' ', zero_effect),
                            'NaN', 'NaN', 'NaN',#bias, rmse, false_positive_perc,
                            'NaN','NaN',#mean(De$StdErr), sd(De$Estimate),
                            mean(De$convergedRuns),
                            'NaN','NaN', #fpp_lower, fpp_upper,
                            num_seedsets,
                            length(De$Estimate), totalRuns,
                            'NaN', # % in CI
                            sep=' & ')
                        cat('\\\\\n')
                        
                        next
                    }
                    
                    totalRuns <- unique(De$totalRuns)
                    stopifnot(length(totalRuns) == 1)
        

                    if (effect == "R_Attribute1"  && substr(error_analysis, 1, 7) == "statnet") {
                        # nodefactor in statnet appears to be not quite the same as R (Activity) in PNet, adjust by subtracting 0.5 from estimate (log-odds) for the statnet value
                        De$Estimate <- De$Estimate - 0.5
                    }

                    rmse <- sqrt(mean((De$Estimate - true_parameters[i])^2))
                    bias <- mean(De$Estimate - true_parameters[i])

                    # count number of times the CI (zSigma std errors) includes true value
                    num_in_ci <- sum((De$Estimate < true_parameters[i] & De$Estimate + zSigma*De$StdErr >= true_parameters[i]) | (De$Estimate >= true_parameters[i] & De$Estimate - zSigma*De$StdErr <= true_parameters[i]))
                    perc_in_ci <- 100 * num_in_ci / length(De$Estimate)

                    # for purely inference (sign and significance, not actual
                    # value of estimate) compute False Postivie rate, as the
                    # number of times the estimate CI does not include zero
                    if (zero_effect == effect) {
                        stopifnot(true_parameters[i] == 0.0)
                        false_positive_count <- sum(
                            (De$Estimate < 0 & De$Estimate + zSigma*De$StdErr < 0) |
                            (De$Estimate > 0 & De$Estimate - zSigma*De$StdErr > 0) )
                        false_positive_perc <- 100 * false_positive_count / length(De$Estimate)
                        confint <- scoreci(false_positive_count, length(De$Estimate), 0.95)
                        fpp_lower <- confint$conf.int[1] * 100
                        fpp_upper <- confint$conf.int[2] * 100

                    }
                    else {
                        false_positive_perc <- NA
                        fpp_lower <- NA
                        fpp_upper <- NA
                    }


                    ## print(results_filename)#xxx
                    ## print(effect)#xxx
                    ## print(zero_effect)#xxx
                    ## print(De$N)#xxx
                    ## print(num_seedsets)

                    cat(network_N, attribute_descr, error_analysis, fixeddensity, 
                        gsub('_', ' ', effect), sub('_', ' ', zero_effect),
                        bias, rmse, false_positive_perc,
                        mean(De$StdErr), sd(De$Estimate), mean(De$convergedRuns),
                        fpp_lower, fpp_upper, num_seedsets,
                        length(De$Estimate), totalRuns, perc_in_ci,
                        sep=' & ')
                    cat('\\\\\n')
		}


}

