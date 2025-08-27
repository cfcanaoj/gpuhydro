// C++17 / MPI / OpenACC
// STL (std::vector / iostream) not used
// Comments are written in English

#include <mpi.h>
#include <openacc.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace mpiconfig_mod {

// ===== Global state =====
static const int mreq = 300;
static MPI_Status stat_[mreq];
static MPI_Request req_[mreq];

int ierr = 0, myid_w = 0, nprocs_w = 0;
MPI_Comm mpi_comm_hyd = MPI_COMM_NULL;
int myid_hyd = 0, nprocs_hyd = 0;
MPI_Comm comm3d = MPI_COMM_NULL;
int myid = 0, nprocs = 0;

int periodic[3] = {1, 1, 1};  // broadcasted later
int ntiles[3]   = {1, 2, 2};
int coords[3]   = {0, 0, 0};
int reorder     = 0;

int n1m=-1, n1p=-1, n2m=-1, n2p=-1, n3m=-1, n3p=-1;
int nreq = 0, nsub = 0;
int gpuid = -1, ngpus = 0;

#pragma acc declare create(myid_w)

// ---- Initialize MPI world, split communicator, create 3D Cartesian topology ----
inline void InitializeMPI() {
  MPI_Init(nullptr, nullptr);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs_w);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid_w);

  if (myid_w == 0) {
    std::printf("MPI process= %d\n", nprocs_w);
    std::printf("decomposition= %d %d %d\n", ntiles[0], ntiles[1], ntiles[2]);
  }
  MPI_Bcast(ntiles, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(periodic, 3, MPI_INT, 0, MPI_COMM_WORLD);

  int np_hyd = ntiles[0] * ntiles[1] * ntiles[2];
  int color  = (np_hyd > 0) ? (myid_w / np_hyd) : 0;
  int key    = myid_w;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &mpi_comm_hyd);
  MPI_Comm_size(mpi_comm_hyd, &nprocs_hyd);
  MPI_Comm_rank(mpi_comm_hyd, &myid_hyd);

  int dims[3] = {ntiles[0], ntiles[1], ntiles[2]};
  int periods[3] = {periodic[0], periodic[1], periodic[2]};
  MPI_Cart_create(mpi_comm_hyd, 3, dims, periods, reorder, &comm3d);
  MPI_Comm_rank(comm3d, &myid);
  MPI_Comm_size(comm3d, &nprocs);

  MPI_Cart_shift(comm3d, 0, 1, &n1m, &n1p);
  MPI_Cart_shift(comm3d, 1, 1, &n2m, &n2p);
  MPI_Cart_shift(comm3d, 2, 1, &n3m, &n3p);
  int c3[3] = {0,0,0};
  MPI_Cart_coords(comm3d, myid, 3, c3);
  coords[0]=c3[0]; coords[1]=c3[1]; coords[2]=c3[2];

  // Select GPU device
  ngpus = acc_get_num_devices(acc_device_nvidia);
  if (myid_w == 0) std::printf("num of GPUs = %d\n", ngpus);
  gpuid = (ngpus > 0) ? (myid_w % ngpus) : -1;
  if (gpuid >= 0) acc_set_device_num(gpuid, acc_device_nvidia);

  #pragma acc update device(myid_w)
}

// ---- Finalize MPI ----
inline void FinalizeMPI() {
  MPI_Finalize();
}

// ---- MINLOC/MAXLOC reductions with CUDA-aware MPI ----
// We use MPI_DOUBLE_INT instead of Fortran's MPI_2DOUBLE_PRECISION
inline void MPIminfind(const double bufinp[2], double bufout[2]) {
  struct Pair { double v; int loc; } in_, out_;
  in_.v = bufinp[0]; in_.loc = (int)bufinp[1];
  #pragma acc host_data use_device(bufinp, bufout)
  {
    MPI_Allreduce(&in_, &out_, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm3d);
  }
  bufout[0] = out_.v; bufout[1] = (double)out_.loc;
}

inline void MPImaxfind(const double bufinp[2], double bufout[2]) {
  struct Pair { double v; int loc; } in_, out_;
  in_.v = bufinp[0]; in_.loc = (int)bufinp[1];
  #pragma acc host_data use_device(bufinp, bufout)
  {
    MPI_Allreduce(&in_, &out_, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm3d);
  }
  bufout[0] = out_.v; bufout[1] = (double)out_.loc;
}

} // namespace mpiconfig_mod

// ==================================================
// DATA IO (MPI-IO)
// ==================================================
namespace mpiiomod {

int ntotal[3] = {0,0,0};   // global dimensions
int npart[3]  = {0,0,0};   // local partition sizes
int nvars = 0, nvarg = 0;  // number of variables

// Buffers (Fortran column-major order assumed)
double* gridX = nullptr; 
double* gridY = nullptr;
double* gridZ = nullptr;
double* data3D = nullptr;  // nvars * npart[0]*npart[1]*npart[2]

static MPI_Datatype SAG1D = MPI_DATATYPE_NULL;
static MPI_Datatype SAG2D = MPI_DATATYPE_NULL;
static MPI_Datatype SAG3D = MPI_DATATYPE_NULL;
static MPI_Datatype SAD3D = MPI_DATATYPE_NULL;

static int is_inited = 0;

static const char id[] = "DT";
static const char datadir[] = "bindata/";

// ---- Helper: local length with edge for last rank ----
inline int local_len_with_edge(int dim) {
  int len = npart[dim];
  if (mpiconfig_mod::coords[dim] == mpiconfig_mod::ntiles[dim]-1) len += 1;
  return len;
}

// ---- malloc wrapper with abort ----
inline void* xmalloc(size_t nbytes) {
  void* p = std::malloc(nbytes);
  if (!p) {
    std::fprintf(stderr, "malloc failed (%zu bytes)\n", nbytes);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return p;
}

// ---- allocate grid buffers ----
inline void allocate_grid_buffers() {
  int li = local_len_with_edge(0);
  int lj = local_len_with_edge(1);
  int lk = local_len_with_edge(2);
  gridX = (double*)xmalloc(sizeof(double) * (size_t)nvarg * (size_t)li);
  gridY = (double*)xmalloc(sizeof(double) * (size_t)nvarg * (size_t)lj);
  gridZ = (double*)xmalloc(sizeof(double) * (size_t)nvarg * (size_t)lk);
}

// ---- allocate data3D buffer ----
inline void allocate_data3D() {
  size_t sz = (size_t)nvars * (size_t)npart[0] * (size_t)npart[1] * (size_t)npart[2];
  data3D = (double*)xmalloc(sizeof(double) * sz);
}

// ---- Output binary files using MPI-IO ----
inline void MPIOutputBindary(int timeid) {
  using namespace mpiconfig_mod;
  int Asize[4], Ssize[4], Start[4];
  MPI_Offset idisp = 0;

  // --- Grid 1D ---
  if (!is_inited) {
    Asize[0]=nvarg; Ssize[0]=nvarg; Start[0]=0;
    Asize[1]=ntotal[0]+1; Ssize[1]=local_len_with_edge(0); Start[1]=npart[0]*coords[0];

    MPI_Type_create_subarray(2, Asize, Ssize, Start, MPI_ORDER_FORTRAN, MPI_DOUBLE, &SAG1D);
    MPI_Type_commit(&SAG1D);

    char fpath[256];
    std::snprintf(fpath, sizeof(fpath), "%sg1d%s", datadir, id);

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAG1D, (char*)"NATIVE", MPI_INFO_NULL);
    MPI_Offset count = (MPI_Offset)Ssize[1] * (MPI_Offset)nvarg;
    MPI_File_write_all(fh, gridX, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  // --- Grid 2D ---
  if (!is_inited) {
    Asize[0]=nvarg; Ssize[0]=nvarg; Start[0]=0;
    Asize[1]=ntotal[1]+1; Ssize[1]=local_len_with_edge(1); Start[1]=npart[1]*coords[1];

    MPI_Type_create_subarray(2, Asize, Ssize, Start, MPI_ORDER_FORTRAN, MPI_DOUBLE, &SAG2D);
    MPI_Type_commit(&SAG2D);

    char fpath[256];
    std::snprintf(fpath, sizeof(fpath), "%sg2d%s", datadir, id);

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAG2D, (char*)"NATIVE", MPI_INFO_NULL);
    MPI_Offset count = (MPI_Offset)Ssize[1] * (MPI_Offset)nvarg;
    MPI_File_write_all(fh, gridY, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  // --- Grid 3D ---
  if (!is_inited) {
    Asize[0]=nvarg; Ssize[0]=nvarg; Start[0]=0;
    Asize[1]=ntotal[2]+1; Ssize[1]=local_len_with_edge(2); Start[1]=npart[2]*coords[2];

    MPI_Type_create_subarray(2, Asize, Ssize, Start, MPI_ORDER_FORTRAN, MPI_DOUBLE, &SAG3D);
    MPI_Type_commit(&SAG3D);

    char fpath[256];
    std::snprintf(fpath, sizeof(fpath), "%sg3d%s", datadir, id);

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAG3D, (char*)"NATIVE", MPI_INFO_NULL);
    MPI_Offset count = (MPI_Offset)Ssize[1] * (MPI_Offset)nvarg;
    MPI_File_write_all(fh, gridZ, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  // --- data3D datatype ---
  if (!is_inited) {
    Asize[0]=nvars; Ssize[0]=nvars; Start[0]=0;
    Asize[1]=ntotal[0]; Ssize[1]=npart[0]; Start[1]=npart[0]*coords[0];
    Asize[2]=ntotal[1]; Ssize[2]=npart[1]; Start[2]=npart[1]*coords[1];
    Asize[3]=ntotal[2]; Ssize[3]=npart[2]; Start[3]=npart[2]*coords[2];

    MPI_Type_create_subarray(4, Asize, Ssize, Start, MPI_ORDER_FORTRAN, MPI_DOUBLE, &SAD3D);
    MPI_Type_commit(&SAD3D);
  }

  // --- write data3D ---
  {
    char filename[64];
    std::snprintf(filename, sizeof(filename), "d3d%s.%05d", id, timeid);
    char fpath[256];
    std::snprintf(fpath, sizeof(fpath), "%s%s", datadir, filename);

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAD3D, (char*)"NATIVE", MPI_INFO_NULL);

    MPI_Offset count =
      (MPI_Offset)nvars *
      (MPI_Offset)npart[0] *
      (MPI_Offset)npart[1] *
      (MPI_Offset)npart[2];

    MPI_File_write_all(fh, data3D, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  is_inited = 1;
}

} // namespace mpiiomod

// ---- Example usage (commented out) ----
// int main() {
//   mpiconfig_mod::InitializeMPI();
//   // Set mpiiomod::nvarg, mpiiomod::nvars, mpiiomod::ntotal[], mpiiomod::npart[]
//   // mpiiomod::allocate_grid_buffers();
//   // mpiiomod::allocate_data3D();
//   // Fill arrays
//   // mpiiomod::MPIOutputBindary(0);
//   // std::free(mpiiomod::gridX); std::free(mpiiomod::gridY);
//   // std::free(mpiiomod::gridZ); std::free(mpiiomod::data3D);
//   mpiconfig_mod::FinalizeMPI();
//   return 0;
// }
