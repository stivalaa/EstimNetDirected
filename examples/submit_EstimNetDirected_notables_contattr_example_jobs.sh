#!/bin/bash

# Submit job array slurm script to do batches of EstimNetDirected runs
# wait for them to finish and then submit job to process output with R scripts

jobid=$(sbatch --parsable run_EstimNetDirected_notables_example_contattr_tnt_slurm_script.sh)

echo submitted job ${jobid}

sbatch --dependency=afterok:${jobid} run_R_scripts_example_contattr_tnt_slurm_script.sh

