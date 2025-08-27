
/**
 * @file boundary.cpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <mpi.h>
#include <omp.h>

#include "hydro_arrays.hpp"
#include "resolution.hpp"
#include "hydro.hpp"
#include "mpi_routines.hpp"

using namespace hydro_arrays_mod;

namespace boundary_mod {
#pragma omp declare target
  Boundary3D<double> Bs,Br;
#pragma omp end declare target 
};

using namespace hydflux_mod;

void AllocateBoundaryVariables(Boundary3D<double>& Bs,Boundary3D<double>& Br){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  Bs.allocate(nprim ,ngh,ktot,jtot,itot);

#pragma omp target update to (Bs.n1, Bs.n2, Bs.n3, Bs.ng, Bs.nv)
#pragma omp target enter data map (alloc: Bs.Xs_data[0:Bs.size1], Bs.Xe_data[0: Bs.size1])
#pragma omp target enter data map (alloc: Bs.Ys_data[0:Bs.size2], Bs.Ye_data[0: Bs.size2])
#pragma omp target enter data map (alloc: Bs.Zs_data[0:Bs.size3], Bs.Ze_data[0: Bs.size3])
  
  Br.allocate(nprim ,ngh,ktot,jtot,itot);

#pragma omp target update to (Br.n1, Br.n2, Br.n3, Br.ng, Br.nv)
#pragma omp target enter data map (alloc: Br.Xs_data[0:Br.size1], Br.Xe_data[0: Br.size1])
#pragma omp target enter data map (alloc: Br.Ys_data[0:Br.size2], Br.Ye_data[0: Br.size2])
#pragma omp target enter data map (alloc: Br.Zs_data[0:Br.size3], Br.Ze_data[0: Br.size3])

}


void DeallocateBoundaryVariables(Boundary3D<double>& Bs,Boundary3D<double>& Br){
  using namespace resolution_mod;
  using namespace hydflux_mod;


#pragma omp target exit data map (delete: Bs.Xs_data[0:Bs.size1], Bs.Xe_data[0: Bs.size1])
#pragma omp target exit data map (delete: Bs.Ys_data[0:Bs.size2], Bs.Ye_data[0: Bs.size2])
#pragma omp target exit data map (delete: Bs.Zs_data[0:Bs.size3], Bs.Ze_data[0: Bs.size3])
  
#pragma omp target exit data map (delete: Bs.Xs_data[0:Br.size1], Br.Xe_data[0: Br.size1])
#pragma omp target exit data map (delete: Bs.Ys_data[0:Br.size2], Br.Ye_data[0: Br.size2])
#pragma omp target exit data map (delete: Bs.Zs_data[0:Br.size3], Br.Ze_data[0: Br.size3])
}


void SendRecvBoundary(const Boundary3D<double>& Bs,Boundary3D<double>& Br){
  using namespace mpiconfig_mod;
  using namespace resolution_mod;
  const int dev = omp_get_default_device();
  void* d_Bs_Xs = omp_get_mapped_ptr(Bs.Xs_data, dev);
  void* d_Bs_Xe = omp_get_mapped_ptr(Bs.Xe_data, dev);
  void* d_Br_Xs = omp_get_mapped_ptr(Br.Xs_data, dev);
  void* d_Br_Xe = omp_get_mapped_ptr(Br.Xe_data, dev);
  
  void* d_Bs_Ys = omp_get_mapped_ptr(Bs.Ys_data, dev);
  void* d_Bs_Ye = omp_get_mapped_ptr(Bs.Ye_data, dev);
  void* d_Br_Ys = omp_get_mapped_ptr(Br.Ys_data, dev);
  void* d_Br_Ye = omp_get_mapped_ptr(Br.Ye_data, dev);
  
  void* d_Bs_Zs = omp_get_mapped_ptr(Bs.Zs_data, dev);
  void* d_Bs_Ze = omp_get_mapped_ptr(Bs.Ze_data, dev);
  void* d_Br_Zs = omp_get_mapped_ptr(Br.Zs_data, dev);
  void* d_Br_Ze = omp_get_mapped_ptr(Br.Ze_data, dev);

  int rc;
  const int dir1=0,dir2=1,dir3=2;
  int nreq = 0;
  
  if(ntiles[dir1] == 1){    
  // |     |Bs.Xe   Bs.Xs|     |
  // |Br.Xs|             |Br.Xs|
#pragma omp target teams distribute parallel for collapse(4)
    for (int n=0; n<nprim; n++)
      for (int k=ks; k<=ke; k++)
	for (int j=js; j<=je;j++)
	  for (int i=0; i<ngh; i++) {
	    Br.Xs(n,k,j,i) = Bs.Xs(n,k,j,i);
	    Br.Xe(n,k,j,i) = Bs.Xe(n,k,j,i);
	  }
  } else {
  // |     |Bs.Xe   Bs.Xs|     |
  // |Br.Xs|             |Br.Xs|
    //#pragma omp target data use_device_addr(Br.Xs_data,Br.Xe_data,Bs.Xs_data,Bs.Xe_data)
#pragma omp target data use_device_addr(d_Bs_Xs,d_Bs_Xe,d_Br_Xs,d_Br_Xe)
    {
    nreq = nreq + 1;
    rc = MPI_Irecv(Br.Xs_data,Br.size1, MPI_DOUBLE, n1m, 1100, comm3d, &req[nreq]);
       
    nreq = nreq + 1;
    rc = MPI_Isend(Bs.Xe_data,Bs.size1, MPI_DOUBLE, n1m, 1200, comm3d, &req[nreq]);

    nreq = nreq + 1;
    rc = MPI_Irecv(Br.Xs_data,Br.size1, MPI_DOUBLE, n1p, 1200, comm3d, &req[nreq]);

    nreq = nreq + 1;
    rc = MPI_Isend(Bs.Xs_data,Bs.size1, MPI_DOUBLE, n1p, 1100, comm3d, &req[nreq]);
    }
   }


  if(ntiles[dir2] == 1){    
  // |     |Bs.Ye   Bs.Ys|     |
  // |Br.Ys|             |Br.Ys|
#pragma omp target teams distribute parallel for collapse(4)
    for (int n=0; n<nprim; n++)
      for (int k=ks; k<=ke; k++)
	for (int j=0; j<ngh;j++)
	  for (int i=is; i<=ie; i++) {
	    Br.Ys(n,k,j,i) = Bs.Ys(n,k,j,i);
	    Br.Ye(n,k,j,i) = Bs.Ye(n,k,j,i);
	  }
  } else {
  // |     |Bs.Ye   Bs.Ys|     |
  // |Br.Ys|             |Br.Ys|
    //#pragma omp target data use_device_addr(Br.Ys_data,Br.Ye_data,Bs.Ys_data,Bs.Ye_data)
#pragma omp target data use_device_addr(d_Bs_Ys,d_Bs_Ye,d_Br_Ys,d_Br_Ye)
    {
    nreq = nreq + 1;
    rc = MPI_Irecv(Br.Ys_data,Br.size2, MPI_DOUBLE, n2m, 2100, comm3d, &req[nreq]);
       
    nreq = nreq + 1;
    rc = MPI_Isend(Bs.Ye_data,Bs.size2, MPI_DOUBLE, n2m, 2200, comm3d, &req[nreq]);

    nreq = nreq + 1;
    rc = MPI_Irecv(Br.Ys_data,Br.size2, MPI_DOUBLE, n2p, 2200, comm3d, &req[nreq]);

    nreq = nreq + 1;
    rc = MPI_Isend(Bs.Ys_data,Bs.size2, MPI_DOUBLE, n2p, 2100, comm3d, &req[nreq]);
    }
   }


  if(ntiles[dir3] == 1){    
  // |     |Bs.Ze   Bs.Zs|     |
  // |Br.Zs|             |Br.Zs|
#pragma omp target teams distribute parallel for collapse(4)
    for (int n=0; n<nprim; n++)
      for (int k=0; k<ngh; k++)
	for (int j=js; j<=je;j++)
	  for (int i=is; i<=ie; i++) {
	    Br.Zs(n,k,j,i) = Bs.Zs(n,k,j,i);
	    Br.Ze(n,k,j,i) = Bs.Ze(n,k,j,i);
	  }
  } else {
  // |     |Bs.Ze   Bs.Zs|     |
  // |Br.Zs|             |Br.Zs|
    //#pragma omp target data use_device_addr(Br.Zs_data,Br.Ze_data,Bs.Zs_data,Bs.Ze_data)
#pragma omp target data use_device_addr(d_Bs_Zs,d_Bs_Ze,d_Br_Zs,d_Br_Ze)
    {
    nreq = nreq + 1;
    rc = MPI_Irecv(Br.Zs_data,Br.size3, MPI_DOUBLE, n3m, 3100, comm3d, &req[nreq]);
       
    nreq = nreq + 1;
    rc = MPI_Isend(Bs.Ze_data,Bs.size3, MPI_DOUBLE, n3m, 3200, comm3d, &req[nreq]);

    nreq = nreq + 1;
    rc = MPI_Irecv(Br.Zs_data,Br.size3, MPI_DOUBLE, n3p, 3200, comm3d, &req[nreq]);

    nreq = nreq + 1;
    rc = MPI_Isend(Bs.Zs_data,Bs.size3, MPI_DOUBLE, n3p, 3100, comm3d, &req[nreq]);
    }
   }
	
  if(nreq != 0) MPI_Waitall ( nreq, req, stat);
  nreq = 0;  
}


void SetBoundaryCondition(Array4D<double>& P,Boundary3D<double>& Bs,Boundary3D<double>& Br) {
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace hydflux_mod;
  using namespace mpiconfig_mod;

  // |     |Bs.Xe   Bs.Xs|     |
  // |Br.Xs|             |Br.Xs|          
  
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=0; i<ngh; i++) {
	  Bs.Xs(n,k,j,i) = P(n,k,j,ie-ngh+1+i);
	  Bs.Xe(n,k,j,i) = P(n,k,j,is      +i);
	  //Xs.data[((n*ktot+k)*jtot+j)*ngh+i]=P.data[((n*ktot+k)*jtot+j)*itot+ie-ngh+1+i];
	  //Xe.data[((n*ktot+k)*jtot+j)*ngh+i]=P.data[((n*ktot+k)*jtot+j)*itot+is      +i];
  }
  
  // y-direction
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=0; j<ngh;j++)
	for (int i=is; i<=ie; i++) {
	  Bs.Ys(n,k,j,i) = P(n,k,je-ngh+1+j,i);
	  Bs.Ye(n,k,j,i) = P(n,k,js      +j,i);
	  //Ys.data[((n*ktot+k)*ngh+j)*itot+i]=P.data[((n*ktot+k)*jtot+je-ngh+1+j)*itot+i];
	  //Ye.data[((n*ktot+k)*ngh+j)*itot+i]=P.data[((n*ktot+k)*jtot+js      +j)*itot+i];
  }
  
  // z-direction
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  Bs.Zs(n,k,j,i) = P(n,ke-ngh+1+k,j,i);
	  Bs.Ze(n,k,j,i) = P(n,ks      +k,j,i);
	  //Zs.data[((n*ngh+k)*jtot+j)*itot+i]=P.data[((n*ktot+ke-ngh+1+k)*jtot+j)*itot+i];
	  //Ze.data[((n*ngh+k)*jtot+j)*itot+i]=P.data[((n*ktot+ks      +k)*jtot+j)*itot+i];
  }

  SendRecvBoundary(Bs,Br);
  
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=0; i<ngh; i++) {
	  P(n,k,j,is-ngh+i) = Br.Xs(n,k,j,i);
	  P(n,k,j,ie+1  +i) = Br.Xe(n,k,j,i);
	  //P.data[((n*ktot+k)*jtot+j)*itot+is-ngh+i] = Xs.data[((n*ktot+k)*jtot+j)*ngh+i];
	  //P.data[((n*ktot+k)*jtot+j)*itot+ie+1  +i] = Xe.data[((n*ktot+k)*jtot+j)*ngh+i];
  }

#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=0; j<ngh;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,k,js-ngh+j,i) = Br.Ys(n,k,j,i);
	  P(n,k,je+1  +j,i) = Br.Ye(n,k,j,i);
	  //P.data[((n*ktot+k)*jtot+js-ngh+j)*itot+i] = Ys.data[((n*ktot+k)*ngh+j)*itot+i];
	  //P.data[((n*ktot+k)*jtot+je+1  +j)*itot+i] = Ye.data[((n*ktot+k)*ngh+j)*itot+i];
  }

#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,ks-ngh+k,j,i) = Br.Zs(n,k,j,i);
	  P(n,ke+1  +k,j,i) = Br.Ze(n,k,j,i);
	  //P.data[((n*ktot+ks-ngh+k)*jtot+j)*itot+i] = Zs.data[((n*ngh+k)*jtot+j)*itot+i];
	  //P.data[((n*ktot+ke+1  +k)*jtot+j)*itot+i] = Ze.data[((n*ngh+k)*jtot+j)*itot+i];
  }
};
