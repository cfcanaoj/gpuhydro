#! /bin/bash
#SBATCH --partition=dgx-full
#SBATCH --gres=gpu:1

# usage sbatch sj_g00.sh
# other useful commands
# sinfo
# squeue

module load nvhpc
nsys profile -o khprof ./kh.x > log.dat
