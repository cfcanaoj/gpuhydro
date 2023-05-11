#PBS -N DT
#PBS -j oe
#PBS -q bulk-cfca
#PBS -l nodes=1
#PBS -m b

# Go to this job's working directory
cd  $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
date  >& log.$PBS_JOBID
aprun -n 8 ./Simulation.x  >> log.$PBS_JOBID
date  >> log.$PBS_JOBID
