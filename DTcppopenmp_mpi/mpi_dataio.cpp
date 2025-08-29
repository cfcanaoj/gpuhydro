
#include<string>
#include "mpi_config.hpp"
#include "mpi_dataio.hpp"

// ==================================================
// DATA IO (MPI-IO)
// ==================================================
namespace mpi_dataio_mod {

  int ntotal[3] = {0,0,0};
  int npart[3]  = {0,0,0};
  int nvars = 0, nvarg = 0;

  GridArray<double> gridXout;
  GridArray<double> gridYout;
  GridArray<double> gridZout;

  FieldArray<double> Fieldout;
  
  static MPI_Datatype SAG1D = MPI_DATATYPE_NULL;
  static MPI_Datatype SAG2D = MPI_DATATYPE_NULL;
  static MPI_Datatype SAG3D = MPI_DATATYPE_NULL;
  static MPI_Datatype SAD3D = MPI_DATATYPE_NULL;
  
  static bool is_inited = false;

  char id[10];
  char datadir[10];

// ---- MPI-IO 書き出し本体（Fortran: MPIOutputBindary）----
  void MPIOutputBindary(int timeid) {
  using namespace mpi_config_mod;

  int Asize[4], Ssize[4], Start[4];
  MPI_Offset idisp = 0;

  // ============ Grid 1D ============
  if (!is_inited) {
    Asize[0] = nvarg; Ssize[0] = nvarg; Start[0] = 0;
    Asize[1] = ntotal[dir1] + 1;
    Ssize[1] = ntotal[dir1];
    if(coords[dir1]==ntiles[dir1]-1) Ssize[1] += 1;
    Start[1] = npart[dir1] * coords[dir1];

    MPI_Type_create_subarray(2, Asize, Ssize, Start, MPI_ORDER_C,MPI_DOUBLE, &SAG1D);
    MPI_Type_commit(&SAG1D);

    std::string fpath = std::string(datadir) + "g1d" + id;
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAG1D, const_cast<char*>("NATIVE"), MPI_INFO_NULL);

    const MPI_Offset count = static_cast<MPI_Offset>(Ssize[1]) * nvarg;
    MPI_File_write_all(fh, gridXout.data, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  // ============ Grid 2D ============
  if (!is_inited) {
    Asize[0] = nvarg; Ssize[0] = nvarg; Start[0] = 0;
    Asize[1] = ntotal[dir2] + 1;
    Ssize[1] = ntotal[dir2];
    if(coords[dir2]==ntiles[dir2]-1) Ssize[1] += 1;
    Start[1] = npart[1] * coords[1];

    MPI_Type_create_subarray(2, Asize, Ssize, Start, MPI_ORDER_C,MPI_DOUBLE, &SAG2D);
    MPI_Type_commit(&SAG2D);

    std::string fpath = std::string(datadir) + "g2d" + id;
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAG2D, const_cast<char*>("NATIVE"), MPI_INFO_NULL);

    const MPI_Offset count = static_cast<MPI_Offset>(Ssize[1]) * nvarg;
    MPI_File_write_all(fh, gridYout.data, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  // ============ Grid 3D ============
  if (!is_inited) {
    Asize[0] = nvarg; Ssize[0] = nvarg; Start[0] = 0;
    Asize[1] = ntotal[dir3] + 1;
    Ssize[1] = ntotal[dir3];
    if(coords[dir3]==ntiles[dir3]-1) Ssize[1] += 1;
    Start[1] = npart[dir3] * coords[dir3];

    MPI_Type_create_subarray(2, Asize, Ssize, Start, MPI_ORDER_FORTRAN,MPI_DOUBLE, &SAG3D);
    MPI_Type_commit(&SAG3D);

    std::string fpath = std::string(datadir) + "g3d" + id;
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAG3D, const_cast<char*>("NATIVE"), MPI_INFO_NULL);

    const MPI_Offset count = static_cast<MPI_Offset>(Ssize[1]) * nvarg;
    MPI_File_write_all(fh, gridZout.data, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  // ============ Data 3D タイプ準備 ============
  if (!is_inited) {
    Asize[0] = nvars; Ssize[0] = nvars; Start[0] = 0;
    Asize[1] = ntotal[0]; Ssize[1] = npart[0]; Start[1] = npart[0] * coords[0];
    Asize[2] = ntotal[1]; Ssize[2] = npart[1]; Start[2] = npart[1] * coords[1];
    Asize[3] = ntotal[2]; Ssize[3] = npart[2]; Start[3] = npart[2] * coords[2];

    MPI_Type_create_subarray(4, Asize, Ssize, Start, MPI_ORDER_FORTRAN,
                             MPI_DOUBLE, &SAD3D);
    MPI_Type_commit(&SAD3D);
  }

  // ============ Data 書き出し ============
  {
    char filename[64];
    std::snprintf(filename, sizeof(filename), "d3d%s.%05d", id, timeid);
    std::string fpath = std::string(datadir) + filename;

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, idisp, MPI_DOUBLE, SAD3D, const_cast<char*>("NATIVE"), MPI_INFO_NULL);

    const MPI_Offset count =
      static_cast<MPI_Offset>(nvars) *
      static_cast<MPI_Offset>(npart[0]) *
      static_cast<MPI_Offset>(npart[1]) *
      static_cast<MPI_Offset>(npart[2]);

    MPI_File_write_all(fh, Fieldout.data, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
  }

  is_inited = true;
}

} // namespace mpiiomod
