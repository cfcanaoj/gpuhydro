# 3D magneto-hydrodynamic deacaying turbulence

## How to copy the source code
After you login the server, `g00.cfca.nao.ac.jp`, follow the instruction.

    cd /gwork0/<username>
    git clone git@github.com:cfcanaoj/gpuhydro gpuhydro
    cd gpuhydro/DTf90openaccc
    

## How to run

### compile 
To run the code, you need to compile 'Simulation.f90' in GPU server.
    
    make Simulation.x
    
Then `Simulation.x`is made in this directory.

### run
Let's run the code.
    
    sbatch sj_g00.sh
    
The simulation data is saved in `bindata/`.

### analysis
To analyze the data, let us make `Analysis.x` in GPU server..
    
    make Analysis.x
    
Now you have many time-snapshots of data. To count it, use a script.
    
    CountBindata.sh
   
See the file, `cat control.dat`. You can know the number of files.
Then preparation is done. Run the analyis.
    
    sbatch sj_g00_ana.sh
    
The output is saved in `output/`.
### 2D plots and animation.
If you need 2D snapshots. cp `bindata` and output` in GPU server to 'HYD2D' in analyis server.
    
    cp -r ???/bindata .
    cp -r ???/output .
    make 2Dsnaps
   
Using `output/vor*.dat`, image files are made and save as `figures/vor*.png`.
To make movie from the files. Type as follows.

    make movie
   
The movie files in saved in `movie/anivor`.

### spectrum
If you need 2D snapshots. cp `bindata` and output` in GPU server to `HYD2D` in analyis server.
To obtain the spectrum
   
    cp -r ???/bindata .
    cp -r ???/output .
    make spectrum
      
## initial condition
stream function is given.

$$ \psi \propto \sin(2\pi k x)\sin(2\pi k y).$$

Initial velocity and magnetic fields are set from $\psi$.

$$ v_x = v_0 \partial_y \psi, v_y = - v_0 \partial_x \psi, $$

$$ B_x = B_0 \partial_y \psi, B_y = - B_0 \partial_x \psi, $$

## analysis

Animation of vorticity and current density are made.

$$ \omega =  \partial_x v_y - \partial_y v_x $$

$$ j =  \partial_x B_y - \partial_y B_x $$

The magnetic potential,$a$, is also calculated.

$$ \vec{B} =  \nabla a $$

The definition of cross helicity is as follows.

$$ h^c = v\cdot B $$

The spectrum of variable, $X$ ,is calculated as follows. In the analysis, $X=E$, $\omega^2$, $a$ and $h^c$.

First 2D Fourier transformation is given.

$$\hat{X}_{{\rm 2D},c}(k_x,k_y) = \iint dx dy X \cos(2\pi (k_x x+k_y y))$$

$$\hat{X}_{{\rm 2D},s}(k_x,k_y) = \iint dx dy X\sin(2\pi (k_x x+k_y y))$$

Then 1D Fourier transformation is given.

$$\hat{X}_{\rm 1D} dk = \sqrt{ \hat{X}_{{\rm 2D},c}^2+\hat{X}_{{\rm 2D},s}^2} dk_x dk_y $$



