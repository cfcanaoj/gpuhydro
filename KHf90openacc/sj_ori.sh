#! /bin/bash
#SBATCH --partition=ga80-1gpu
#SBATCH --gres=gpu:1

# usage sbatch sj_ori.sh
# other useful commands
# sinfo
# squeue

module load nvhpc

# original program
./kh_ori.x > log.dat
