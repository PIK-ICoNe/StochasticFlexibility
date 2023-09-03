#!/bin/bash

#SBATCH --time=00-24:00:00
#SBATCH --qos=medium
#SBATCH --mem=16000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --array=1-9%9

#SBATCH --output=output/jobarray-%A_%a.out

echo "SLURM TASK ID: $SLURM_ARRAY_TASK_ID"

module load julia

julia ./code/StochasticFlexibility/paper_experiments/conv_simulate.jl new_parameters
