#!/bin/bash

#SBATCH --time=00-36:00:00
#SBATCH --qos=medium
#SBATCH --mem=24000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --array=1-9%9

#SBATCH --output=output/convergence-%A_%a.out

echo "SLURM TASK ID: $SLURM_ARRAY_TASK_ID"

module load julia

julia ./code/StochasticFlexibility/paper_experiments/conv_simulate.jl improved_constraint
