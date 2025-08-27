# gpuhydro
To foster GPU computing, we show sample codes of hydrodynamic simulations in this repositry.

|Problem|Language|parallel|Name|
----|----|----|----
|2D Kelvin-Helmholtz Instability|Fortran|OpenMP(CPU)|[KHf90openmp](KHf90openmp)|
|2D Kelvin-Helmholtz Instability|Fortran|MPI, OpenMP(CPU)|[KHf90openmp_mpi](KHf90openmp_mpi)|
|2D Kelvin-Helmholtz Instability|Fortran|OpenACC|[KHf90openacc](KHf90openacc)|
|2D Kelvin-Helmholtz Instability|Fortran|MPI, OpenACC|[KHf90openacc_mpi](KHf90openacc_mpi)|
|3D upwind advection|C++|OpenMP(GPU)|[ADcppopenmp](ADcppopenmp)|
|3D upwind advection|C++|OpenACC|[ADcppopenacc](ADcppopenacc)|
|3D magneto-hydrodynamic deacaying turbulence|Fortran|MPI, OpenMP(CPU)|[DTf90openmp_mpi](DTf90openmp_mpi)|
|3D magneto-hydrodynamic deacaying turbulence|Fortran|OpenACC|[DTf90openacc](DTf90openacc)|
|3D magneto-hydrodynamic deacaying turbulence|Fortran|MPI, OpenACC|[DTf90openacc_mpi](DTf90openacc_mpi)|
|3D magneto-hydrodynamic deacaying turbulence|C++|OpenMP(GPU)|[DTcppopenmp](DTcppopenmp)|


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
