#!/usr/bin/Rscript
##
## File:    plot_example_fourcycles_dependence_figures.R
## Author:  Alex Stivala
## Created: June 2024
##
## Make figures illustrating dependence structure for node-centered
## four-cycles statistic.
##
## Usage: Rscript plot_example_fourcycles_dependence_figures.R
##
## Writes output files:
##
##    karate_11_1.pdf
##    karate_34_9.pdf
##    southern_women_5_23.pdf
##    southern_women_6_21.pdf
##
## in cwd (WARNING: overwrites).
##
## The CYPATH directory containing cypath executable and Perl scripts
## (transgrh.pl etc.) must be in the PATH as this R script uses
## system2() to call them in order to find (not just count)
## four-cycles. Have to edit the transgrh.pl script to add a
## semicolon on the end of line 38 as it does not work otherwise.
##
## For CYPATH see:
##
##    http://research.nii.ac.jp/~uno/code/cypath.html
##    http://research.nii.ac.jp/~uno/code/cypath11.zip
## 
##    Uno, T., & Satoh, H. (2014, October). An efficient algorithm for
##    enumerating chordless cycles and chordless paths. In International
##    Conference on Discovery Science (pp. 313-324). Springer, Cham.
##
##


library(igraph)
library(RColorBrewer)
### install.packages("remotes")
##remotes::install_github("schochastics/networkdata")
library(networkdata)

## read in R source file from directory where this script is located
##http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep=.Platform$file.sep))
}
source_local('plot_fourcycles_dependence_figure.R')


###
### main
###

print(date())
print(sessionInfo())

## one-mode: Zachary Karate Club
plot_fourcycles_dependence_figure(karate,11,1,'karate_11_1.pdf')
plot_fourcycles_dependence_figure(karate,34,9,'karate_34_9.pdf')

## two-mode (bipartite): Southern Women
plot_fourcycles_dependence_figure(southern_women,5,23,'southern_women_5_23.pdf')
plot_fourcycles_dependence_figure(southern_women,6,21,'southern_women_6_21.pdf')
