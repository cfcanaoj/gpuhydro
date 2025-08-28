#!/usr/bin/env bash
#SBATCH --partition=dgx-full
#SBATCH --nodes=1
#SBATCH --gpu-bind=closest
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --time=00:05:00
#SBATCH -o ./out%j.log
#SBATCH -e ./err%j.log


module purge
module load nvhpc/25.7

#export LIBOMPTARGET_DEBUG=1
#export LIBOMPTARGET_INFO=4

cd ${SLURM_SUBMIT_DIR}
mpiexec -n 2 ./main.x
