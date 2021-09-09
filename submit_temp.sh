#!/bin/bash

#SBATCH --partition=dufek
#SBATCH --account=dufeklab
#SBATCH --job-name=ME_F_R
#SBATCH --output=ME.out
#SBATCH --error=ME.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=300000M

# This post-processing script assumes binary data files






module load matlab/R2019b
matlab -nosplash -nodisplay -singleCompThread -r "temperaturePP" 
###"compPP('FULL_COMP','binary')"





