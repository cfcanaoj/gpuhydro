#ifndef MPI_ROUTINES_HPP
#define MPI_ROUTINES_HPP

#include <mpi.h>
#include <cstddef> // size_t

// -----------------------------
// Namespace: mpimod
// -----------------------------
namespace mpiconfig_mod {

 inline constexpr int mreq = 300;
 inline MPI_Status stat[mreq];
 inline MPI_Request req[mreq];
  
// ---- Public MPI state (read-mostly) ----
extern MPI_Comm mpi_comm_hyd;   // split communicator
extern MPI_Comm comm3d;         // 3D Cartesian communicator

extern int ierr;
extern int myid_w, nprocs_w;    // WORLD rank/size
extern int myid,   nprocs;      // comm3d rank/size

extern int periodic[3];         // periodicity flags (0/1)
extern int ntiles[3];           // process grid dims (x,y,z)
extern int coords[3];           // my coords in the grid (x,y,z)
extern int reorder;             // MPI_Cart_create reorder flag

// neighbor ranks (prev/next) along each dimension
extern int n1m, n1p, n2m, n2p, n3m, n3p;

// GPU info (OpenACC)
extern int gpuid, ngpus;


}

void InitializeMPI();  // init world, split, create Cartesian, select GPU

#endif
