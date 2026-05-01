#!/bin/bash

#SBATCH --partition=katira
#SBATCH -J S4_nb4_af005
#SBATCH --time=14-0:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=6
#SBATCH --mem=58G
#SBATCH --output=%j.output
#SBATCH --error=%j.error

#Load the conda environment:
module load matlab/2023b

# run the MATLAB code
matlab -nodisplay -nosplash -nodesktop -r "run('sim_nb4_af005_N1C00038.m');exit;" &
sleep 300
matlab -nodisplay -nosplash -nodesktop -r "run('sim_nb4_af005_N1C0034.m');exit;" &
sleep 300
matlab -nodisplay -nosplash -nodesktop -r "run('sim_nb4_af005_N25C0047.m');exit;" &
sleep 300
matlab -nodisplay -nosplash -nodesktop -r "run('sim_nb4_af005_N25C1.m');exit;" &
sleep 300
matlab -nodisplay -nosplash -nodesktop -r "run('sim_nb4_af005_N5C003.m');exit;" &
sleep 300
matlab -nodisplay -nosplash -nodesktop -r "run('sim_nb4_af005_N5C03.m');exit;" &
wait
