#!/bin/bash

#SBATCH --partition=dufek     ### Partition
#SBATCH --account=dufeklab
#SBATCH --job-name=si_dev2  ### Job Name
#SBATCH --output=si_dev2.out   ### file in which to store job stdout
#SBATCH --error=si_dev2.err    ### file in which to store job stderr
#SBATCH --time=36:00:00      ### WallTime
#SBATCH --nodes=1           ### Number of Nodes
#SBATCH --ntasks-per-node=40 ### Number of tasks (MPI processes)

####SBATCH --mem-per-cpu=376000M


#./comp1.exe > out1
#wait
#./comp2.exe > out2
#wait
#./comp3.exe > out3
#wait
#./comp4.exe > out4
#wait
#./comp5.exe > out5
#wait
#./comp6.exe > out6
#wait
#./comp11.exe > out11
#wait
echo "ending thermodynamics, starting dynamics"
timeout 20h ./melting.exe
wait
tag=$( tail -n 1 TIMESTEPS )
echo "at end of cycle 1, timesteps = ,$tag"
wait
#while [ $tag -lt 1 ]
#do
#timeout 1h ./melting.exe
#wait
#tag=$( tail -n 1 TIMESTEPS )
#echo "while cycle, timesteps = ,$tag"
done
##wait
##echo "moving to transfer files"
##./transfer.exe

