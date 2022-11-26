# Kelvin-Helmholtz Instability Test

## How to copy the source code
After you login the server, `xc.cfca.nao.ac.jp`, follow the instruction.

    cd /work/<username>
    git clone git@github.com:cfcanaoj/gpuhydro gpuhydro
    cd gpuhydro/KHf90openmp
    

## How to run
Before running job, pleaes check the category of you in `pbs_xc.sh`.

	module switch PrgEnv-cray PrgEnv-intel
	make
	qsub pbs_xc.sh

## How to see the results
Let us move to analysis server.

    ssh an10@cfca.nao.ac.jp
    cp /xc-work/<username>/gpuhydro/KHf90openaccc .
    cd KHf90openaccc
    module load gnuplot
    gnuplot dn2dx.plt
    display figures/dnx00100.png

## Description of the problem

https://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html

## Numerical setup

https://ui.adsabs.harvard.edu/abs/2012ApJS..201...18M/abstract
