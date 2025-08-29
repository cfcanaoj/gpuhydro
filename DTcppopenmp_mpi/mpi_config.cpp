
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <omp.h>

#include "mpi_config.hpp"

// =============================
// Namespace: mpimod (definitions)
// =============================
namespace mpi_config_mod {

  // Public state
  int ierr = 0;
  int myid_w = 0, nprocs_w = 0;
  
  MPI_Comm mpi_comm_hyd = MPI_COMM_NULL;
  int myid_hyd = 0, nprocs_hyd = 0;
  
  MPI_Comm comm3d = MPI_COMM_NULL;
  int myid = 0, nprocs = 0;
  
  int periodic[3] = {0, 0, 0};   // broadcasted as MPI_INT
  int ntiles[3]   = {0, 0, 0};
  int coords[3]   = {0, 0, 0};
  int reorder     = 0;
  int n1m=-1, n1p=-1, n2m=-1, n2p=-1, n3m=-1, n3p=-1;

  int gpuid = -1, ngpus = 0;

  // ---- Initialize MPI world, split communicator, create 3D Cartesian topology ---


void InitializeMPI() {
  MPI_Init(nullptr, nullptr);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs_w);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid_w);
  
  if (myid_w == 0) {
    std::printf("MPI process= %d\n", nprocs_w);
    std::printf("decomposition= %d %d %d\n", ntiles[0], ntiles[1], ntiles[2]);
  }
  MPI_Bcast(  ntiles, ndim, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(periodic, ndim, MPI_INT, 0, MPI_COMM_WORLD);
  
  // Split communicator into groups of size ntiles[0]*ntiles[1]*ntiles[2]
  int np_hyd = ntiles[0] * ntiles[1] * ntiles[2];
  int color  = (np_hyd > 0) ? (myid_w / np_hyd) : 0;
  int key    = myid_w;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &mpi_comm_hyd);
  MPI_Comm_size(mpi_comm_hyd, &nprocs_hyd);
  MPI_Comm_rank(mpi_comm_hyd, &myid_hyd);
    
  // Create 3D Cartesian topology
  MPI_Cart_create(mpi_comm_hyd, ndim, ntiles, periodic, reorder, &comm3d);
  MPI_Comm_rank(comm3d, &myid);
  MPI_Comm_size(comm3d, &nprocs);

  // Neighbors and my coordinates
  MPI_Cart_shift(comm3d, dir1, 1, &n1m, &n1p);
  MPI_Cart_shift(comm3d, dir2, 1, &n2m, &n2p);
  MPI_Cart_shift(comm3d, dir3, 1, &n3m, &n3p);
  MPI_Cart_coords(comm3d, myid, 3,coords);

  
  // Select GPU device (OpenMP)
  ngpus = omp_get_num_devices();
  if (myid_w == 0) std::printf("num of GPUs = %d\n", ngpus);
  gpuid = (ngpus > 0) ? (myid_w % ngpus) : -1;
  if (gpuid >= 0) omp_set_default_device(gpuid);

  printf("init myid_w,gpuid=(%i,%i)\n",myid_w,gpuid);  
    // If device-side use of myid_w is needed later, this updates it once.
#pragma acc update device(myid_w)
}

void MPIminfind(const double& vin,const int& locin, double& vout, int& locout){
  struct Pair { double v; int loc; } in_, out_;
  in_.v   =   vin;
  in_.loc = locin;
  MPI_Allreduce(&in_,& out_, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm3d);
    vout  = out_.v;
  locout  = out_.loc;
}

void MPImaxfind(const double& vin,const int& locin, double& vout, int& locout) {
  struct Pair { double v; int loc; } in_, out_;
  in_.v   =   vin;
  in_.loc = locin;
  MPI_Allreduce(&in_,&out_, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm3d);
    vout  = out_.v;
  locout  = out_.loc;
}
};
