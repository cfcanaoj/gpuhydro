#!/usr/bin/env bash
#SBATCH --partition=dgx-full
#SBATCH --gres=gpu:1
#SBATCH --time=00:05:00

module purge
module load nvhpc/25.7
#module load gnuplot
cd ${SLURM_SUBMIT_DIR}
./main.x
