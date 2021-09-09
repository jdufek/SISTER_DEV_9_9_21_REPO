#!/bin/bash

#SBATCH --partition=dufek     ### Partition/queue name specific to our group
#SBATCH --job-name=mat_sc  ### Job Name -- can make this specific to your program
#SBATCH --output=mat.out   ### file in which to store job stdout, edit for a your case
#SBATCH --error=mat.err    ### file in which to store job stderr, edit for your case
#SBATCH --time=48:00:00      ### WallTime (maximum running time)
#SBATCH --nodes=1           ### Number of Nodes
#SBATCH --ntasks-per-node=1 ### Number of tasks -- this is set for a single core executable

#####SBATCH --mem-per-cpu=376000M

module load matlab
matlab -nosplash -nodisplay -singleCompThread -r ldm_plots   ### Put your executable name here
