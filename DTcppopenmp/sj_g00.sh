#!/usr/bin/env bash
#SBATCH --partition=ga80-1gpu
#SBATCH --gres=gpu:1
#SBATCH --time=00:05:00
#SBATCH -o ./out%j.log
#SBATCH -e ./err%j.log


module purge
module load nvhpc/25.7

export LIBOMPTARGET_DEBUG=1
export LIBOMPTARGET_INFO=4

cd ${SLURM_SUBMIT_DIR}
./main.x
