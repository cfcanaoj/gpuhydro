#! /bin/bash
#SBATCH --partition=ga80-1gpu
#SBATCH --gres=gpu:1

# usage sbatch sj_g00.sh
# other useful commands
# sinfo
# squeue

module load nvhpc
./Simulation.x > log.dat
