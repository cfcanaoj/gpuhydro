#PBS -N KH
#PBS -j oe
#PBS -q bulk-cfca
#PBS -l nodes=1
#PBS -m b

# Go to this job's working directory
cd  $PBS_O_WORKDIR
export OMP_NUM_THREADS=2

date  >& log.$PBS_JOBID
time aprun -cc none -n 16 ./kh.x  >> log.$PBS_JOBID
date  >> log.$PBS_JOBID
