#!/usr/bin/Rscript
##
## File:    plot_asymmetric_BipartiteFourCyclesNodePower_simulation_results_uniquenodes.R
## Author:  Alex Stivala
## Created: April 2025
##
## Plot box plots of statistics of simulated bipartite networks
## (generated with
## run_SimulateERGM_bipartite_Example_A12000_B4000_attrs_tnt_slurm_script.sh
## script), but this time plotting the number of unique A nodes (left axis
## red) and unique B nodes (right axis blue) on the same plot.
## These counts are read from the output of the
## run_make_unique_bipartite_nodes_count_table.sh
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

## Note second axis can only be a linear transform of first axis
## assuming here that number of B nodes is half number of A nodes
## as it is in the asymmetric example (100 A  and 50 B)

stopifnot(dat$numnodesB == dat$numnodesA/2)

p <- ggplot(data = dat[dat$theta_BipartiteFourCyclesNodePowerAxB != "neg.neg",],
            aes(x = theta_BipartiteFourCyclesNodePowerAxB)) +
  geom_boxplot(aes(y = uniqueFourCyclesA), color='red') +
  geom_boxplot(aes(y = uniqueFourCyclesB), color='blue') +
  scale_y_continuous(name = 'Unique A nodes in four-cycles',
                     sec.axis = sec_axis(trans = ~./2,
                                         name='Unique B nodes in four-cycles')) +
  theme_classic() + theme(axis.title.y = element_text(colour = 'red'),
                          axis.title.y.right = element_text(colour = 'blue'))

postscript("bipartite_asymmetric_fourcyclesnodepower_simulation_uniquenodes_boxplots.eps")
print(p)
dev.off()


