#!/usr/bin/Rscript
##
## File:    plot_asymmetric_BipartiteFourCyclesNodePower_simulation_results.R
## Author:  Alex Stivala
## Created: September 2024
##
## Plot box plots of statistics of simulated bipartite networks
## (generated with
## run_SimulateERGM_asymmetric_BipartiteFourCyclesNodePower_examples.sh
## script).
##
## Output file in cwd (WARNING: overwrite if exists):
##
##   bipartite_asymmetric_fourcyclesnodepower_simulation_boxplots.eps
##

library(ggplot2)
library(grid)
library(gridExtra)

statfiles <- list('stats_sim_bipartite_asymmetric_FourCyclesNodePower_Anegative_Bnegative.txt',
                  'stats_sim_bipartite_asymmetric_FourCyclesNodePower_Anegative_Bpositive.txt',
                  'stats_sim_bipartite_asymmetric_FourCyclesNodePower_Apositive_Bnegative.txt',
                  'stats_sim_bipartite_asymmetric_FourCyclesNodePower_Apositive_Bpositive.txt',
                  'stats_sim_bipartite_asymmetric_FourCyclesNodePower_Apositive_Bzero.txt',
                  'stats_sim_bipartite_asymmetric_FourCyclesNodePower_Azero_Bpositive.txt',
                  'stats_sim_bipartite_asymmetric_FourCyclesNodePower_Azero_Bzero.txt')
stats_dfs <- lapply(statfiles, function(x) read.table(x, header=T))
stats_dfs[[1]]$theta_BipartiteFourCyclesNodePowerA <- "neg"
stats_dfs[[2]]$theta_BipartiteFourCyclesNodePowerA <- "neg"
stats_dfs[[3]]$theta_BipartiteFourCyclesNodePowerA <- "pos"
stats_dfs[[4]]$theta_BipartiteFourCyclesNodePowerA <- "pos"
stats_dfs[[5]]$theta_BipartiteFourCyclesNodePowerA <- "pos"
stats_dfs[[6]]$theta_BipartiteFourCyclesNodePowerA <- "zero"
stats_dfs[[7]]$theta_BipartiteFourCyclesNodePowerA <- "zero"
stats_dfs[[1]]$theta_BipartiteFourCyclesNodePowerB <- "neg"
stats_dfs[[2]]$theta_BipartiteFourCyclesNodePowerB <- "pos"
stats_dfs[[3]]$theta_BipartiteFourCyclesNodePowerB <- "neg"
stats_dfs[[4]]$theta_BipartiteFourCyclesNodePowerB <- "pos"
stats_dfs[[5]]$theta_BipartiteFourCyclesNodePowerB <- "zero"
stats_dfs[[6]]$theta_BipartiteFourCyclesNodePowerB <- "pos"
stats_dfs[[7]]$theta_BipartiteFourCyclesNodePowerB <- "zero"

stats <- rbind(stats_dfs[[1]], stats_dfs[[2]], stats_dfs[[3]],
               stats_dfs[[4]], stats_dfs[[5]], stats_dfs[[6]],
               stats_dfs[[7]])

stats$theta_BipartiteFourCyclesNodePowerAxB <-  interaction(stats$theta_BipartiteFourCyclesNodePowerA, stats$theta_BipartiteFourCyclesNodePowerB)


plotlist <- list()

p <- ggplot(data = stats,
            aes(x = theta_BipartiteFourCyclesNodePowerAxB,
                y = BipartiteFourCyclesNodePowerA.5.)) +
  theme_classic() +
  ylab("BipartiteFourCyclesNodePowerA") +
  geom_boxplot()

plotlist <- c(plotlist, list(p))

p <- ggplot(data = stats,
            aes(x = theta_BipartiteFourCyclesNodePowerAxB,
                y = BipartiteFourCyclesNodePowerB.5.)) +
  theme_classic() +
  ylab("BipartiteFourCyclesNodePowerB") +
  geom_boxplot()

plotlist <- c(plotlist, list(p))


p <- ggplot(data = stats,
            aes(x = theta_BipartiteFourCyclesNodePowerAxB,
                y = FourCycles)) +
  theme_classic() +
  geom_boxplot()

plotlist <- c(plotlist, list(p))


postscript("bipartite_asymmetric_fourcyclesnodepower_simulation_boxplots.eps")
do.call(grid.arrange, plotlist)
dev.off()


