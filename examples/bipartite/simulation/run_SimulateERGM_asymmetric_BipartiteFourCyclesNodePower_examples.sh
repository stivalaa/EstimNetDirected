#!/bin/sh

echo -n "started at: "; date
uname -a

ROOT=../../..

for configfile in config_sim_bipartite_asymmetric_FourCyclesNodePower_Anegative_Bnegative.txt config_sim_bipartite_asymmetric_FourCyclesNodePower_Anegative_Bpositive.txt config_sim_bipartite_asymmetric_FourCyclesNodePower_Apositive_Bnegative.txt config_sim_bipartite_asymmetric_FourCyclesNodePower_Apositive_Bpositive.txt config_sim_bipartite_asymmetric_FourCyclesNodePower_Apositive_Bzero.txt config_sim_bipartite_asymmetric_FourCyclesNodePower_Azero_Bpositive.txt config_sim_bipartite_asymmetric_FourCyclesNodePower_Azero_Bzero.txt
do
  echo ========== ${configfile} ==========
  time ${ROOT}/src/SimulateERGM  ${configfile}
  # Last awk removes leading and trailing whitespace on filename, see:
  # https://unix.stackexchange.com/questions/102008/how-do-i-trim-leading-and-trailing-whitespace-from-each-line-of-some-output
  statsfile=`grep statsFile ${configfile} | cut -d= -f2 | awk '{$1=$1};1'`
  time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R ${statsfile}
done

Rscript plot_asymmetric_BipartiteFourCyclesNodePower_simulation_results.R
Rscript plot_asymmetric_BipartiteFourCyclesNodePower_simulation_graph_visualizations.R

# Also plot 'goodness-of-fit' plots where 'observed' network is the
# single simulated network used in graph visualizations

for prefix in simulation_sim_bipartite_asymmetric_FourCyclesNodePower_Anegative_Bnegative simulation_sim_bipartite_asymmetric_FourCyclesNodePower_Anegative_Bpositive simulation_sim_bipartite_asymmetric_FourCyclesNodePower_Apositive_Bnegative simulation_sim_bipartite_asymmetric_FourCyclesNodePower_Apositive_Bpositive simulation_sim_bipartite_asymmetric_FourCyclesNodePower_Apositive_Bzero simulation_sim_bipartite_asymmetric_FourCyclesNodePower_Azero_Bpositive simulation_sim_bipartite_asymmetric_FourCyclesNodePower_Azero_Bzero
do
  echo ========== ${prefix} ==========
  Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R -t -y 16 ${prefix}_9900000.net ${prefix}
done

echo -n "ended at: "; date

