#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --qos=priority
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --array=1-2%2
#SBATCH --output=output/jobarray-%A_%a.out


# 1. make sure all cores we get are on one node
#    jobs of this type may not span nodes.


# 2. we're a 1-task job


# 3. this is how many CPUs (threads) our task will need
#    (up to 16 (Haswell) or 32 (Broadwell))


module load julia

julia ./code/StochasticFlexibility/paper_experiments/baseline_simulate.jl new_constraint
