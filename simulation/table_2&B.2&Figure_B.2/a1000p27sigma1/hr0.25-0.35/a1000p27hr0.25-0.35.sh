#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --array=1-10
#SBATCH --job-name="a1000p27hr0.25-0.35"
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%a.out

# Some points on the above
# job-name   : this just specifies what job is called, here I have set it so the name of the job is the function being run
#              in the actual script
# array      : currently set 1-1, but can set as x-y (x,y integers x<y), this then runs the same job a number of times with indexes
#              x,x+1,x+2,...,y and the index takes on the value $SLURM_ARRAY_TASK_ID
# mail-user  : this option sends an email to the address written once job is done

# Load required modules and set directory
module load StdEnv/2020 r/4.1.0

# This runs the script; here instead of manually entering the 5 arguments, the first 4 are based on what we already defined in the sbatch
# commands above here
# $SLURM_JOB_NAME         = function being run in Rscript
# $SLURM_ARRAY_TASK_ID    = job number
# $SLURM_CPUS_PER_TASK    = number of cores to use
# $SLURM_TASK_COUNT       = 1 as only running one task
# CEDARSCRATCH            = name of the compute server

# Define parameters within the script
index=$SLURM_ARRAY_TASK_ID          # Batch index (1-10)

# Define other parameters (these values are fixed for all jobs)
model="a"                           # Model type ("a" or "b")
sample_size=1000                      # Sample size
dimension=27                        # Dimension
candidate="40:200"                   # Candidate vector (in R format)
bandwidth="0.25 0.3 0.35"  # Bandwidth vector
GCV="FALSE"                          # GCV parameter (TRUE or FALSE)
simulation=40                       # Simulation for each batch (e.g., 40 for batch 1-10 means 400 simulations in total)


# Print the input values for debugging purposes
echo "Batch index: $index"
echo "Model: $model"
echo "Sample size: $sample_size"
echo "Dimension: $dimension"
echo "Candidate vector: $candidate"
echo "Bandwidth: $bandwidth"
echo "GCV: $GCV"  # Print the GCV value
echo "Simulation: $simulation"  # Print the simulation value

Rscript ../../simulation_table2.R $index $model $sample_size $dimension "$candidate" "$bandwidth" "$GCV" $simulation
