#!/bin/bash
#SBATCH --job-name=KH
#SBATCH --partition=M-test-cfca
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH -o slmlog%J.out
#SBATCH --hint=nomultithread

cd ${SLURM_SUBMIT_DIR}

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

date >& log.$SLURM_JOB_ID
#time srun ./kh.x >> log.$SLURM_JOB_ID
time srun vtune -collect hotspots -r=vtout$SLURM_JOB_ID  ./kh.x >> log.$SLURM_JOB_ID 
date >> log.$SLURM_JOB_ID

