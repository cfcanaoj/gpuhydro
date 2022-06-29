# Kelvin-Helmholtz Instability Test

## How to copy the source code

    cd /gwork0/<username>
    git clone git@github.com:cfcanaoj/gpuhydro gpuhydro
    cd gpuhydro/KHf90openaccc
    

## How to run

	module load nvhpc
	make
	sbatch sj_g00.sh

## How to see the results
Let us move to analysis server.

    ssh an10@cfca.nao.ac.jp
    cd /gwork0/<username>/gpuhydro/KHf90openaccc
    module load gnuplot
    gnuplot dn2dx.plt
    display figures/dnx00100.png

## Description of the problem

https://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html

## Numerical setup

https://ui.adsabs.harvard.edu/abs/2012ApJS..201...18M/abstract



