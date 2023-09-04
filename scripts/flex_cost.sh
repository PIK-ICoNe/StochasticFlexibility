#!/bin/bash

#SBATCH --time=00-08:00:00
#SBATCH --qos=medium
#SBATCH --mem=16000
#SBATCH --output=output/jobarray-%A_%a.out
#SBATCH --array=1-18%16

#SBATCH --nodes=1

#SBATCH --ntasks=1

echo "SLURM TASK ID: $SLURM_ARRAY_TASK_ID"

module load julia

srun julia ./code/StochasticFlexibility/paper_experiments/flex_cost.jl new_constraint 3 3



