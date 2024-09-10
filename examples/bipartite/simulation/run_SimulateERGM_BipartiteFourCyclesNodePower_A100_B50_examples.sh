#!/bin/sh

echo -n "started at: "; date
uname -a

ROOT=../../..

for configfile in config_sim_bipartite_A100_B50_FourCyclesNodePower_Anegative_Bnegative.txt config_sim_bipartite_A100_B50_FourCyclesNodePower_Anegative_Bpositive.txt config_sim_bipartite_A100_B50_FourCyclesNodePower_Apositive_Bnegative.txt config_sim_bipartite_A100_B50_FourCyclesNodePower_Apositive_Bpositive.txt config_sim_bipartite_A100_B50_FourCyclesNodePower_Apositive_Bzero.txt config_sim_bipartite_A100_B50_FourCyclesNodePower_Azero_Bpositive.txt config_sim_bipartite_A100_B50_FourCyclesNodePower_Azero_Bzero.txt
do
  echo ========== ${configfile} ==========
  time ${ROOT}/src/SimulateERGM  ${configfile}
  # Last awk removes leading and trailing whitespace on filename, see:
  # https://unix.stackexchange.com/questions/102008/how-do-i-trim-leading-and-trailing-whitespace-from-each-line-of-some-output
  statsfile=`grep statsFile ${configfile} | cut -d= -f2 | awk '{$1=$1};1'`
  time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R ${statsfile}
done

Rscript plot_BipartiteFourCyclesNodePower_simulation_results.R

echo -n "ended at: "; date

