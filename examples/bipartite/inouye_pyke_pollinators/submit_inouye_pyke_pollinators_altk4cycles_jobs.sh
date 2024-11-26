#!/bin/bash

# Submit estimation job and subsequet jobs to do convergence and GoF
# checks, using slurm job dependencies to only run after the jobs they 
# depend on finish successfully.

estim_jobid=$(sbatch --parsable run_EstimNetDirected_inouye_pyke_pollinators_altk4cycles_slurm_script.sh)
echo ${estim_jobid}
# "pseudo-gof" convergence check and simulation (for degeneracy and gof check)
# can run simultaneously, after estimation finished
sbatch --dependency=afterok:${estim_jobid} run_plotEstimNetDirectedSimFit_inouye_pyke_pollinators_altk4cycles_R_slurm_script.sh 
sim_jobid=$(sbatch --parsable --dependency=afterok:${estim_jobid} run_SimulatERGM_gof_inouye_pyke_pollinators_altk4cycles_slurm_script.sh )
echo ${sim_jobid}
# GoF plot depends on output of simulation
sbatch --dependency=afterok:${sim_jobid} run_plot_GofSimFit_inouye_pyke_pollinators_altk4cycles_R_slurm_script.sh 
