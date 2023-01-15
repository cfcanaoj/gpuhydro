#! /bin/bash
#SBATCH --partition=ga80-1gpu
#SBATCH --gres=gpu:1

# usage sbatch sj_g00.sh
# other useful commands
# sinfo
# squeue

module load nvhpc

# original program
./kh_ori.x > log.dat

# practice program
#./kh_pra.x > log.dat
