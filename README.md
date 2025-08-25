# gpuhydro
To foster GPU computing, we show sample codes of hydrodynamic simulations in this repositry.

|Problem|Language|directive|Name|
----|----|----|----
|3D upwind advection|C++|GPU, OpenMP|[ADcppopenmp](ADcppopenmp)|
|3D upwind advection|C++|GPU, OpenACC|[ADcppopenacc](ADcppopenacc)|
|3D magneto-hydrodynamic deacaying turbulence|Fortran|GPU, OpenACC|[DTf90openacc](DTf90openacc)|
|3D magneto-hydrodynamic deacaying turbulence|C++|GPU, OpenMP|[DTcppopenmp](DTcppopenmp)|
|2D Kelvin-Helmholtz Instability|Fortran|GPU, OpenACC|[KHf90openacc](KHf90openacc)|
|2D Kelvin-Helmholtz Instability|Fortran|CPU, OpenMP|[KHf90openmp](KHf90openmp)|

## Fortran with OpenACC

### 2D Kelvin-Helmholtz Instability
The code for GPU computaiton is prepared in [KHf90openacc](KHf90openacc). To compare the performance with its CPU version, we have also prepared  openmp version in [KHf90openmp](KHf90openmp).

### 3D magneto-hydrodynamic deacaying turbulence
The code for GPU computaiton is prepared in [DTf90openacc](DTf90openacc).

## C++ with OpenMP

### 3D upwind advection
The code for GPU computaiton is prepared in [ADcppopenmp](ADcppopenmp).

### 3D magneto-hydrodynamic deacaying turbulence
The code for GPU computaiton is prepared in [DTcppopenmp](DTcppopenmp).

# Link
- [CfCA GPU Cluster](https://www.cfca.nao.ac.jp/gpgpu)
