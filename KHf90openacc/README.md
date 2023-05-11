# Kelvin-Helmholtz Instability Test

## How to copy the source code
First you need the setup for ssh-connection to github (see [the instrcution](../HowToUseGithub.md)).

After you login the server, `g00.cfca.nao.ac.jp`, perform the following command.
    
    cd /gwork0/<username>
    git clone -b ws2022 git@github.com:cfcanaoj/gpuhydro gpuhydro
    cd gpuhydro/KHf90openaccc
    

## How to run

	module load nvhpc
	make
	sbatch sj_ori.sh

## How to see the results
Let us move to analysis server.

    ssh <username>@an.cfca.nao.ac.jp
    cd /gwork0/<username>/gpuhydro/KHf90openaccc
    module load gnuplot
    gnuplot dn2dx_ori.plt
    display figures_ori/dnx00070.png

## Description of the problem

https://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html

## Profile
Run the code with profiler.

    sbatch sj_g00prof.sh

In `sj_g00prof.sh' the profiler is called as follows. 
   
   nsys profile -o khprof ./kh_ori.x
   
The profile data is summarized in `khprof.nsys-rep`.

## Numerical setup

https://ui.adsabs.harvard.edu/abs/2012ApJS..201...18M/abstract



