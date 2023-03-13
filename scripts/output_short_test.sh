#!/bin/bash

#SBATCH --time=00-00:05:00
#SBATCH --output=output/%x-%j.out
# The simplest of all SLURM submission scripts.
# Submit with "sbatch ex1.sh".

# We'll run the following command on the first (and in this case, only)
# task in the allocation of resources SLURM gives us.

module load julia

srun julia ./code/StochasticFlexibility/manuscripts/test_flexibility.jl
hostname

# The output will be in a file called 'slurm-<jobid>.out' and will
# look something like:
# 'cs-f14c05b04'
# i.e. the hostname of the server the script ran on.
