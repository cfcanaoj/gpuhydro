#! /bin/bash
#SBATCH --partition=ga80-1gpu
#SBATCH --gres=gpu:1
#SBATCH -o ./out%j.log
#SBATCH -e ./err%j.log

# usage sbatch sj_ori.sh
# other useful commands
# sinfo
# squeue

module load nvhpc

# original program
./kh.x
