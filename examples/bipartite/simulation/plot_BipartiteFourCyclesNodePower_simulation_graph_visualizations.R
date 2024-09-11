#!/usr/bin/Rscript
##
## File:    plotBipartiteFourCyclesNodePower_simulation_graph_visualizations.R
## Author:  Alex Stivala
## Created: September 2024
##
## Plot visualizations of simulated networks
## (generated with
## run_SimulateERGM_bipartite_Example_A12000_B4000_attrs_tnt_slurm_script.sh
## script).
##
## Output file in cwd (WARNING: overwrite if exists):
##
##   bipartite_fourcyclesnodepower_simulation_visualizations.eps
##

library(igraph)

plot_graph <- function(g, name) {
  plot(g, vertex.label = NA,
       vertex.size = 4,
       vertex.color = ifelse(V(g)$type, 'blue', 'red'),
       vertex.shape = ifelse(V(g)$type, 'square', 'circle'),
       layout = layout.kamada.kawai,
       main = name
       )
}

graph_prefixes <- list(
  'simulation_sim_bipartite_FourCyclesNodePower_Anegative_Bnegative',
  'simulation_sim_bipartite_FourCyclesNodePower_Anegative_Bpositive',
  'simulation_sim_bipartite_FourCyclesNodePower_Apositive_Bnegative',
  'simulation_sim_bipartite_FourCyclesNodePower_Apositive_Bpositive',
  'simulation_sim_bipartite_FourCyclesNodePower_Apositive_Bzero',
  'simulation_sim_bipartite_FourCyclesNodePower_Azero_Bpositive',
  'simulation_sim_bipartite_FourCyclesNodePower_Azero_Bzero')

graph_filenames <- lapply(graph_prefixes,
                           function(s) paste(s, '_9900000.net', sep=''))

## lines up with filenames above
graphnames <- list(
  'Anegative_Bnegative',
  'Anegative_Bpositive',
  'Apositive_Bnegative',
  'Apositive_Bpositive',
  'Apositive_Bzero',
  'Azero_Bpositive',
  'Azero_Bzero'
)
stopifnot(length(graphnames) == length(graph_filenames))

postscript("bipartite_fourcyclesnodepower_simulation_visualizations.eps")
par(mfrow = c(2,4),
    mar=c(1, 0, 1, 0)) 
for (i in 1:length(graph_filenames)) {
  filename <- graph_filenames[[i]]
  g <- read.graph(filename, format='pajek')
  print(filename)
  print(graphnames[[i]])
  summary(g)
  plot_graph(g, graphnames[[i]])
}
dev.off()
