# 3D upwind advection

## How to copy the source code
After you login the server, `g00.cfca.nao.ac.jp`, follow the instruction.

    cd /gwork0/<username>
    git clone git@github.com:cfcanaoj/gpuhydro gpuhydro
    cd gpuhydro/ADcppopenmp
    

## How to run

### compile 
To run the code, you need to compile 'Simulation.f90' in GPU server.
    
    make main.x
    
Then `Simulation.x`is made in this directory.

### run
Let's run the code.
    
    sbatch sj_g00.sh
    
The simulation data is saved in `bindata/`.
