#!/bin/bash

#SBATCH --time=00-72:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --array=1-75%16
#SBATCH --partition=largemem

#SBATCH --output=output/jobarray-%A_%a.out

echo "SLURM TASK ID: $SLURM_ARRAY_TASK_ID"

module load julia

julia ./StochasticFlexibility/paper_experiments/conv_simulate.jl run_02_14
