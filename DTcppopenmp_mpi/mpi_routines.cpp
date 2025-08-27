
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <openacc.h>

#include "mpi_routines.hpp"

// =============================
// Namespace: mpimod (definitions)
// =============================
namespace mpiconfig_mod {

  static const int mreq = 300;
  static MPI_Status stat_[mreq];
  static MPI_Request req_[mreq];

  // Public state
  int ierr = 0, myid_w = 0, nprocs_w = 0;
  MPI_Comm mpi_comm_hyd = MPI_COMM_NULL;
  int myid_hyd = 0, nprocs_hyd = 0;
  MPI_Comm comm3d = MPI_COMM_NULL;
  int myid = 0, nprocs = 0;

  int periodic[3] = {1, 1, 1};   // broadcasted as MPI_INT
  int ntiles[3]   = {1, 2, 2};
  int coords[3]   = {0, 0, 0};
  int reorder     = 0;

  int n1m=-1, n1p=-1, n2m=-1, n2p=-1, n3m=-1, n3p=-1;

  int gpuid = -1, ngpus = 0;

  // ---- Initialize MPI world, split communicator, create 3D Cartesian topology ---
};

using namespace mpiconfig_mod;

void InitializeMPI() {
  MPI_Init(nullptr, nullptr);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs_w);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid_w);
  
  if (myid_w == 0) {
    std::printf("MPI process= %d\n", nprocs_w);
    std::printf("decomposition= %d %d %d\n", ntiles[0], ntiles[1], ntiles[2]);
  }
  MPI_Bcast(ntiles, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(periodic, 3, MPI_INT, 0, MPI_COMM_WORLD);
  
  // Split communicator into groups of size ntiles[0]*ntiles[1]*ntiles[2]
  int np_hyd = ntiles[0] * ntiles[1] * ntiles[2];
  int color  = (np_hyd > 0) ? (myid_w / np_hyd) : 0;
  int key    = myid_w;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &mpi_comm_hyd);
  MPI_Comm_size(mpi_comm_hyd, &nprocs_hyd);
  MPI_Comm_rank(mpi_comm_hyd, &myid_hyd);
    
  // Create 3D Cartesian topology
  int dims[3] = {ntiles[0], ntiles[1], ntiles[2]};
  int periods[3] = {periodic[0], periodic[1], periodic[2]};
  MPI_Cart_create(mpi_comm_hyd, 3, dims, periods, reorder, &comm3d);
  MPI_Comm_rank(comm3d, &myid);
  MPI_Comm_size(comm3d, &nprocs);

  // Neighbors and my coordinates
  MPI_Cart_shift(comm3d, 0, 1, &n1m, &n1p);
  MPI_Cart_shift(comm3d, 1, 1, &n2m, &n2p);
  MPI_Cart_shift(comm3d, 2, 1, &n3m, &n3p);
  int c3[3] = {0,0,0};
  MPI_Cart_coords(comm3d, myid, 3, c3);
  coords[0]=c3[0]; coords[1]=c3[1]; coords[2]=c3[2];
  
  // Select GPU device (OpenACC)
  ngpus = acc_get_num_devices(acc_device_nvidia);
  if (myid_w == 0) std::printf("num of GPUs = %d\n", ngpus);
  gpuid = (ngpus > 0) ? (myid_w % ngpus) : -1;
  if (gpuid >= 0) acc_set_device_num(gpuid, acc_device_nvidia);
  
    // If device-side use of myid_w is needed later, this updates it once.
#pragma acc update device(myid_w)
}


