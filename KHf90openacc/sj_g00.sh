#! /bin/bash
#SBATCH --partition=ga80-1gpu
#SBATCH --ntasks=4
#SBATCH --gres=gpu:1

# usage sbatch sj_ori.sh
# other useful commands
# sinfo
# squeue

module load nvhpc

# program
mpiexec -n 4 ./kh.x > log.dat