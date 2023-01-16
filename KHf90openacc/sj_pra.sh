#! /bin/bash
#SBATCH --partition=ga80-1gpu
#SBATCH --gres=gpu:1

# usage sbatch sj_pra.sh
# other useful commands
# sinfo
# squeue

module load nvhpc

# practice program
./kh_pra.x > log.dat
