#!/bin/bash

#SBATCH --qos=medium
#SBATCH --time=00-48:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --array=41-57:4
#SBATCH --partition=largemem

#SBATCH --output=output/jobarray-%A_%a.out

echo "SLURM TASK ID: $SLURM_ARRAY_TASK_ID"

module load julia

julia ./code/StochasticFlexibility/experiments/investment_sensitivity_single_sample.jl
