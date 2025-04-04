#!/usr/bin/Rscript
##
## File:    plot_asymmetric_BipartiteFourCyclesNodePower_simulation_results_uniquenodes.R
## Author:  Alex Stivala
## Created: April 2025
##
## Plot box plots of statistics of simulated bipartite networks
## (generated with
## run_SimulateERGM_bipartite_Example_A12000_B4000_attrs_tnt_slurm_script.sh
## script), but this time plotting the number of unique A nodes
## and unique B nodes on the same plot.
## These counts are read from the output of the
## run_make_asymmetric_unique_bipartite_nodes_count_table.sh
## script.
##
## Output file in cwd (WARNING: overwrite if exists):
##
##   bipartite_asymmetric_fourcyclesnodepower_simulation_uniquenodes_boxplots.eps
##

library(ggplot2)

setEPS()  # postscript() will use EPS settings

dat <- read.table('unique_bipartite_nodes_count_table_bipartite_asymmetric.txt',header=T, stringsAsFactors=TRUE)

dat$theta_BipartiteFourCyclesNodePowerAxB<- interaction(dat$Aparam, dat$Bparam)



p <- ggplot(data = dat[dat$theta_BipartiteFourCyclesNodePowerAxB != "neg.neg",],
            aes(x = theta_BipartiteFourCyclesNodePowerAxB)) +
  geom_boxplot(aes(y = uniqueFourCyclesA), color='red') +
  geom_boxplot(aes(y = uniqueFourCyclesB), color='blue') +
  ylab('Unique A (red) or B (blue) nodes in four-cycles') +
  theme_classic() + 
  xlab("BipartiteFourCyclesNodePower A and B parameter signs") 

postscript("bipartite_asymmetric_fourcyclesnodepower_simulation_uniquenodes_boxplots.eps")
print(p)
dev.off()

