#!/bin/bash

#SBATCH --qos=medium
#SBATCH --time=03-00:00:00
#SBATCH --output=output/%x-%j.out


#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --array=0-1023%32

module load julia

srun julia ./code/StochasticFlexibility/experiments/investment_sensitivity_parallel.jl
