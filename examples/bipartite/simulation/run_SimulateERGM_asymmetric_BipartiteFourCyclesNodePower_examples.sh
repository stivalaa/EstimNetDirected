#!/bin/bash

#SBATCH --job-name="SimulateERGM_asymmetric_BipartiteFourCyclesNodePower"
#SBATCH --ntasks=1
#SBATCH --time=0-02:00:00
#SBATCH --mem-per-cpu=2GB
#SBATCH --output=SimulateERGM_asymmetric_BipartiteFourCyclesNodePower-%j.out
#SBATCH --error=SimulateERGM_asymmetric_BipartiteFourCyclesNodePower-%j.err

echo -n "started at: "; date
uname -a

module load gcc/11.3.0 # needed by r/4.2.1
module load openmpi/4.1.4 # needed by r/4.2.1
module load r/4.2.1

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
  ## no -t as Opsahl cc too slow to compute on (at least) simulation_sim_bipartite_asymmetric_FourCyclesNodePower_Apositive_Bnegative
  ## -y 16 takes too long (not finished after 45 minutes on positive_positive) so does 12 (1.4 hours just for observed), changed to 8 (as even 10 was 421943596)
  Rscript ${ROOT}/scripts/plotEstimNetDirectedSimFit.R  -y 8 ${prefix}_9900000.net ${prefix}
done

echo -n "ended at: "; date

