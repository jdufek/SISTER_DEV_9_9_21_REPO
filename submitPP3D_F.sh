#!/bin/bash

#SBATCH --partition=dufek
#SBATCH --account=dufeklab
#SBATCH --job-name=ldmPP
#SBATCH --output=ldmPP3D.out
#SBATCH --error=ldmPP3D.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=376000M

# This submission script is for the post-processing of Laguna del Maule 
# magmatic simulation outputs unrelated to the forward modeling routines.
# It is tailored for submitting to the U. of Oregon Talapas cluster.

# Parameters for the post-processing routines:

# Timestep filtering: 1 = No; 2 = Yes
filt=1


## WRITE SECONDARY PARAMETER INPUT FILE

# Input filename
INPUTFW=INPUTP

echo "TIMESTEP FILTERING (INTEGER):" > $INPUTFW
echo $filt >> $INPUTFW


## PREPARE FOR POST-PROCESSING

# Make a directory for post-processing products
mkdir POST

# Directory location of scripts
SCRIPTS=/projects/dufeklab/dufek/ldmScripts

# Copy over scripts
PPSCRIPTS=$SCRIPTS/postProcessing

cp $PPSCRIPTS/ldmGeneralPP.m ./
cp $PPSCRIPTS/ldmInputProcessing.m ./

# Load MATLAB and submit post-processing jobs
module load matlab/R2019b


## DENSITY POST-PROCESSING

# Run the post-processing
matlab -nodisplay -batch "ldmGeneralPP('drho'); exit"

# Change permissions and move post-processing products
chmod 775 density3D*
mv density3D*.png density3D*.avi ./POST


## MELT FRACTION POST-PROCESSING

# Run the post-processing
matlab -nodisplay -batch "ldmGeneralPP('f'); exit"

# Change permissions and move post-processing products
chmod 775 mFraction3D*
mv mFraction3D*.png mFraction3D*.avi ./POST

# Breakup the output of the Stdout file
echo " "


## TEMPERATURE POST-PROCESSING

# Run the post-processing
matlab -nodisplay -batch "ldmGeneralPP('temp'); exit"

# Change permissions and move post-processing products
chmod 775 temperature3D*
mv temperature3D*.png temperature3D*.avi ./POST

# Breakup the output of the Stdout file
echo " "


## TEMPERATURE ANOMALY POST-PROCESSING

# Run the post-processing
matlab -nodisplay -batch "ldmGeneralPP('tempanomaly'); exit"

# Change permissions and move post-processing products
chmod 775 tAnomaly3D*
mv tAnomaly3D*.png tAnomaly3D*.avi ./POST
