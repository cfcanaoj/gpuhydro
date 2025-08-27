#PBS -N KH
#PBS -j oe
#PBS -q single
#PBS -l nodes=1:ppn=48
#PBS -m b


module load intel/compiler intel/mpi
# Go to this job's working directory
cd  $PBS_O_WORKDIR
export OMP_NUM_THREADS=12
date  >& log.$PBS_JOBID
mpiexec -n 4 ./kh.x  >> log.$PBS_JOBID
date  >> log.$PBS_JOBID
