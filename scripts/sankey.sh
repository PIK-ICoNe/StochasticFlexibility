#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=00-00:30:00
#SBATCH --output=output/%x-%j.out


# 1. make sure all cores we get are on one node
#    jobs of this type may not span nodes.
#SBATCH --nodes=1

# 2. we're a 1-task job
#SBATCH --ntasks=1

# 3. this is how many CPUs (threads) our task will need
#    (up to 16 (Haswell) or 32 (Broadwell))
#SBATCH --cpus-per-task=1

module load julia

srun julia --threads=1 ./code/StochasticFlexibility/experiments/cost_reduction_sankey.jl



