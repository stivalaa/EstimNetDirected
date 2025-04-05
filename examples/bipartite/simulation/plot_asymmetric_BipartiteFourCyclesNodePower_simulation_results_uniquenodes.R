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

## Change default font size to make it larger so readable when included in
## LaTeX documents and reduced in smaller panels
## https://ggplot2.tidyverse.org/articles/faq-customising.html

## using theme_classic() to get no grey background and no gridlines
## for some journal requirements e.g. J. Complex Networks

theme_set(theme_classic(base_size = 12))

# Also increase x axis labels specifically 
theme_update(axis.text.x = element_text(size = 14))


dat <- read.table('unique_bipartite_nodes_count_table_bipartite_asymmetric.txt',header=T, stringsAsFactors=TRUE)

dat$theta_BipartiteFourCyclesNodePowerAxB<- interaction(dat$Aparam, dat$Bparam)

## make data frame with two copies of the data, one for A nodes and one for
## B nodes with unique A and B node counts respectively so that we can use
## the ggplot2 colour or fill for factor A or B (using colour as some boxplots
## have little variance so just black line on fill, using colour instead
## they are red or blue).

xdat <- dat
xdat$nodeset <- 'A'
ydat <- dat
ydat$nodeset <- 'B'
xdat$unique_count <- dat$uniqueFourCyclesA
ydat$unique_count <- dat$uniqueFourCyclesB
dat2 <- rbind(xdat,ydat)
dat2$nodeset <- factor(dat2$nodeset)



## remove neg.neg case to match other plots, these have no 4-cycles anyway
p <- ggplot(data = dat2[dat2$theta_BipartiteFourCyclesNodePowerAxB != "neg.neg",],
            aes(x = theta_BipartiteFourCyclesNodePowerAxB)) +
  geom_boxplot(aes(y = unique_count, color = nodeset)) +
  ylab('Unique A (red) or B (blue) nodes in four-cycles') +
  xlab("BipartiteFourCyclesNodePower A and B parameter signs") +
  scale_colour_manual(values = c('red','blue'), name = "Node set")

postscript("bipartite_asymmetric_fourcyclesnodepower_simulation_uniquenodes_boxplots.eps", horizontal=FALSE, paper="special", width=9, height=6)
print(p)
dev.off()

