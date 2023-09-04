#!/bin/bash

#SBATCH --time=00-12:00:00
#SBATCH --nodes=1
#SBATCH --qos=medium
#SBATCH --ntasks-per-node=16
#SBATCH --array=1-15%15
#SBATCH --mem=16000

#SBATCH --output=output/jobarray-%A_%a.out

echo "SLURM TASK ID: $SLURM_ARRAY_TASK_ID"

module load julia

julia ./code/StochasticFlexibility/paper_experiments/request_rate_simulate.jl new_constraint
