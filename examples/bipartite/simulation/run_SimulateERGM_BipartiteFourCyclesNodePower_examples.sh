#!/bin/bash

#SBATCH --job-name="SimulateERGM_BipartiteFourCyclesNodePower"
#SBATCH --ntasks=1
#SBATCH --time=0-02:00:00
#SBATCH --mem-per-cpu=1GB
#SBATCH --output=SimulateERGM_BipartiteFourCyclesNodePower-%j.out
#SBATCH --error=SimulateERGM_BipartiteFourCyclesNodePower-%j.err

echo -n "started at: "; date
uname -a

module load gcc/11.3.0 # needed by r/4.2.1
module load openmpi/4.1.4 # needed by r/4.2.1
module load r/4.2.1

ROOT=../../..

for configfile in config_sim_bipartite_FourCyclesNodePower_Anegative_Bnegative.txt config_sim_bipartite_FourCyclesNodePower_Anegative_Bpositive.txt config_sim_bipartite_FourCyclesNodePower_Apositive_Bnegative.txt config_sim_bipartite_FourCyclesNodePower_Apositive_Bpositive.txt config_sim_bipartite_FourCyclesNodePower_Apositive_Bzero.txt config_sim_bipartite_FourCyclesNodePower_Azero_Bpositive.txt config_sim_bipartite_FourCyclesNodePower_Azero_Bzero.txt
do
  echo ========== ${configfile} ==========
  time ${ROOT}/src/SimulateERGM  ${configfile}
  # Last awk removes leading and trailing whitespace on filename, see:
  # https://unix.stackexchange.com/questions/102008/how-do-i-trim-leading-and-trailing-whitespace-from-each-line-of-some-output
  statsfile=`grep statsFile ${configfile} | cut -d= -f2 | awk '{$1=$1};1'`
  time Rscript ${ROOT}/scripts/plotSimulationDiagnostics.R ${statsfile}
done

Rscript plot_BipartiteFourCyclesNodePower_simulation_results.R
Rscript plot_BipartiteFourCyclesNodePower_simulation_graph_visualizations.R

# Also plot 'goodness-of-fit' plots where 'observed' network is the
# single simulated network used in graph visualizations

for prefix in simulation_sim_bipartite_FourCyclesNodePower_Anegative_Bnegative simulation_sim_bipartite_FourCyclesNodePower_Anegative_Bpositive simulation_sim_bipartite_FourCyclesNodePower_Apositive_Bnegative simulation_sim_bipartite_FourCyclesNodePower_Apositive_Bpositive simulation_sim_bipartite_FourCyclesNodePower_Apositive_Bzero simulation_sim_bipartite_FourCyclesNodePower_Azero_Bpositive simulation_sim_bipartite_FourCyclesNodePower_Azero_Bzero
do
  echo ========== ${prefix} ==========
  ## -t and -y 16 too slow on higher density / higher four-cycles graphs
  ## even -y 8 too slow on pos.pos
  Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R  -y 6 ${prefix}_9900000.net ${prefix}
done

echo -n "ended at: "; date

