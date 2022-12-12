#!/bin/bash -x

#PJM -x PJM_LLIO_GFSCACHE=/vol0005
#PJM -N KHtest
#PJM -L "node=1"
#PJM -L "rscgrp=small"
#PJM -L "elapse=1:00:00"
#PJM -S
#PJM --mpi "proc=1"
#PJM --mpi "max-proc-per-node=1"

export PARALLEL=48
export OMP_NUM_THREADS=$PARALLEL
export XOS_MMM_L_HPAGE_TYPE=none

# batch scripts for fugaku
# usage:
# pjsub pj_fk.sh
# pjstat
# pjdel

mpiexec -n 1 ./kh.x

