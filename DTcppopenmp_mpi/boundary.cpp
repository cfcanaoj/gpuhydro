
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

#include "resolution.hpp"
#include "hydro.hpp"
#include "boundary.hpp"

#include "mpi_config.hpp"

using namespace hydflux_mod;

namespace boundary_mod {
#pragma omp declare target
  BoundaryArray<double> Bs,Br;
#pragma omp end declare target 

auto assocb = [&](void* host_ptr, size_t bytes, int dev) {
    void* dptr = omp_get_mapped_ptr(host_ptr, dev);
    if (!dptr) {
        std::fprintf(stderr, "mapped_ptr is NULL for %p\n", host_ptr);
    }
    int r = omp_target_associate_ptr(host_ptr, dptr, bytes, /*device_offset=*/0, dev);
    if (r != 0) {
        std::fprintf(stderr, "omp_target_associate_ptr failed (%d)\n", r);
    }
};

void AllocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br){
  using namespace resolution_mod;

  int dev = omp_get_default_device();
  
  // buffer for send
  Bs.allocate(nprim ,ngh,ktot,jtot,itot);

#pragma omp target enter data map (alloc: Bs.Xs_data[0:Bs.size1], Bs.Xe_data[0: Bs.size1])
#pragma omp target enter data map (alloc: Bs.Ys_data[0:Bs.size2], Bs.Ye_data[0: Bs.size2])
#pragma omp target enter data map (alloc: Bs.Zs_data[0:Bs.size3], Bs.Ze_data[0: Bs.size3])
#pragma omp target update to (Bs.n1, Bs.n2, Bs.n3, Bs.ng, Bs.nv, Bs.size1, Bs.size2, Bs.size3)
  /*
  assocb(Bs.Xs_data, sizeof(double)*Bs.size1,dev);
  assocb(Bs.Xe_data, sizeof(double)*Bs.size1,dev);
  assocb(Bs.Ys_data, sizeof(double)*Bs.size2,dev);
  assocb(Bs.Ye_data, sizeof(double)*Bs.size2,dev);
  assocb(Bs.Zs_data, sizeof(double)*Bs.size3,dev);
  assocb(Bs.Ze_data, sizeof(double)*Bs.size3,dev);
  */
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<ngh; i++) {
	  Bs.Xs(n,k,j,i) = 0.0;
	  Bs.Xe(n,k,j,i) = 0.0;
  }
  
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<ngh; j++)
	for (int i=0; i<itot; i++) {
	  Bs.Ys(n,k,j,i) = 0.0;
	  Bs.Ye(n,k,j,i) = 0.0;
  }
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  Bs.Zs(n,k,j,i) = 0.0;
	  Bs.Ze(n,k,j,i) = 0.0;
  }
#pragma omp target update to (Bs.Xs_data[0:Bs.size1],Bs.Xe_data[0:Bs.size1])
#pragma omp target update to (Bs.Ys_data[0:Bs.size2],Bs.Ye_data[0:Bs.size2])
#pragma omp target update to (Bs.Zs_data[0:Bs.size3],Bs.Ze_data[0:Bs.size3])

  // buffer for receive

  Br.allocate(nprim ,ngh,ktot,jtot,itot);

#pragma omp target enter data map (alloc: Br.Xs_data[0:Br.size1], Br.Xe_data[0: Br.size1])
#pragma omp target enter data map (alloc: Br.Ys_data[0:Br.size2], Br.Ye_data[0: Br.size2])
#pragma omp target enter data map (alloc: Br.Zs_data[0:Br.size3], Br.Ze_data[0: Br.size3])
#pragma omp target update to (Br.n1, Br.n2, Br.n3, Br.ng, Br.nv, Br.size1, Br.size2, Br.size3)
  /*
  assocb(Br.Xs_data, sizeof(double)*Br.size1,dev);
  assocb(Br.Xe_data, sizeof(double)*Br.size1,dev);
  assocb(Br.Ys_data, sizeof(double)*Br.size2,dev);
  assocb(Br.Ye_data, sizeof(double)*Br.size2,dev);
  assocb(Br.Zs_data, sizeof(double)*Br.size3,dev);
  assocb(Br.Ze_data, sizeof(double)*Br.size3,dev);
  */
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<ngh; i++) {
	  Br.Xs(n,k,j,i) = 0.0;
	  Br.Xe(n,k,j,i) = 0.0;
  }
  
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<ngh; j++)
	for (int i=0; i<itot; i++) {
	  Br.Ys(n,k,j,i) = 0.0;
	  Br.Ye(n,k,j,i) = 0.0;
  }
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  Br.Zs(n,k,j,i) = 0.0;
	  Br.Ze(n,k,j,i) = 0.0;
  }

#pragma omp target update to (Br.Xs_data[0:Br.size1],Br.Xe_data[0:Br.size1])
#pragma omp target update to (Br.Ys_data[0:Br.size2],Br.Ye_data[0:Br.size2])
#pragma omp target update to (Br.Zs_data[0:Br.size3],Br.Ze_data[0:Br.size3])


};


void DeallocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br){
  using namespace resolution_mod;
  using namespace hydflux_mod;


#pragma omp target exit data map (delete: Bs.Xs_data[0:Bs.size1], Bs.Xe_data[0: Bs.size1])
#pragma omp target exit data map (delete: Bs.Ys_data[0:Bs.size2], Bs.Ye_data[0: Bs.size2])
#pragma omp target exit data map (delete: Bs.Zs_data[0:Bs.size3], Bs.Ze_data[0: Bs.size3])
  
#pragma omp target exit data map (delete: Br.Xs_data[0:Br.size1], Br.Xe_data[0: Br.size1])
#pragma omp target exit data map (delete: Br.Ys_data[0:Br.size2], Br.Ye_data[0: Br.size2])
#pragma omp target exit data map (delete: Br.Zs_data[0:Br.size3], Br.Ze_data[0: Br.size3])
}


void SendRecvBoundary(const BoundaryArray<double>& Bs,BoundaryArray<double>& Br){
  using namespace mpi_config_mod;
  using namespace resolution_mod;
  int dev = omp_get_default_device();
  //printf("myid gpuid=%i %i",myid_w,dev);
  int rc;
  int nreq = 0;
  
  if(ntiles[dir1] == 1){    
  // |     |Bs.Xe   Bs.Xs|     |
  // |Br.Xs|             |Br.Xe|
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
  // |Br.Xs|             |Br.Xe|
    void* d_Bs_Xs = omp_get_mapped_ptr(Bs.Xs_data, dev);
    void* d_Bs_Xe = omp_get_mapped_ptr(Bs.Xe_data, dev);
    void* d_Br_Xs = omp_get_mapped_ptr(Br.Xs_data, dev);
    void* d_Br_Xe = omp_get_mapped_ptr(Br.Xe_data, dev);
    rc = MPI_Irecv(d_Br_Xs,Br.size1, MPI_DOUBLE, n1m, 1100, comm3d, &req[nreq++]);
    rc = MPI_Isend(d_Bs_Xe,Bs.size1, MPI_DOUBLE, n1m, 1200, comm3d, &req[nreq++]);
    rc = MPI_Irecv(d_Br_Xe,Br.size1, MPI_DOUBLE, n1p, 1200, comm3d, &req[nreq++]);
    rc = MPI_Isend(d_Bs_Xs,Bs.size1, MPI_DOUBLE, n1p, 1100, comm3d, &req[nreq++]);
   
   }


  if(ntiles[dir2] == 1){    
  // |     |Bs.Ye   Bs.Ys|     |
  // |Br.Ys|             |Br.Ye|
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
  // |Br.Ys|             |Br.Ye|
  //               |     |Bs.Ye   Bs.Ys|     |
  //               |Br.Ys|             |Br.Ye|
    void* d_Bs_Ys = omp_get_mapped_ptr(Bs.Ys_data, dev);
    void* d_Bs_Ye = omp_get_mapped_ptr(Bs.Ye_data, dev);
    void* d_Br_Ys = omp_get_mapped_ptr(Br.Ys_data, dev);
    void* d_Br_Ye = omp_get_mapped_ptr(Br.Ye_data, dev);
    
    rc = MPI_Irecv(d_Br_Ys,Br.size2, MPI_DOUBLE, n2m, 2100, comm3d, &req[nreq++]);
    rc = MPI_Isend(d_Bs_Ye,Bs.size2, MPI_DOUBLE, n2m, 2200, comm3d, &req[nreq++]);
    rc = MPI_Irecv(d_Br_Ye,Br.size2, MPI_DOUBLE, n2p, 2200, comm3d, &req[nreq++]);
    rc = MPI_Isend(d_Bs_Ys,Bs.size2, MPI_DOUBLE, n2p, 2100, comm3d, &req[nreq++]);

  }


  if(ntiles[dir3] == 1){    
  // |     |Bs.Ze   Bs.Zs|     |
  // |Br.Zs|             |Br.Ze|
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
  // |Br.Zs|             |Br.Ze|  
    void* d_Bs_Zs = omp_get_mapped_ptr(Bs.Zs_data, dev);
    void* d_Bs_Ze = omp_get_mapped_ptr(Bs.Ze_data, dev);
    void* d_Br_Zs = omp_get_mapped_ptr(Br.Zs_data, dev);
    void* d_Br_Ze = omp_get_mapped_ptr(Br.Ze_data, dev);
    
    rc = MPI_Irecv(d_Br_Zs,Br.size3, MPI_DOUBLE, n3m, 3100, comm3d, &req[nreq++]);
    rc = MPI_Isend(d_Bs_Ze,Bs.size3, MPI_DOUBLE, n3m, 3200, comm3d, &req[nreq++]);
    rc = MPI_Irecv(d_Br_Ze,Br.size3, MPI_DOUBLE, n3p, 3200, comm3d, &req[nreq++]);
    rc = MPI_Isend(d_Bs_Zs,Bs.size3, MPI_DOUBLE, n3p, 3100, comm3d, &req[nreq++]);
   }

    if(nreq != 0) MPI_Waitall ( nreq, req, MPI_STATUSES_IGNORE);
    nreq = 0;	
};


void SetBoundaryCondition(FieldArray<double>& P,BoundaryArray<double>& Bs,BoundaryArray<double>& Br) {
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace hydflux_mod;
  using namespace mpi_config_mod;

  // x-direction
  // |     |Bs.Xe   Bs.Xs|     |
  // |Br.Xs|             |Br.Xe|
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=0; i<ngh; i++) {
	  Bs.Xs(n,k,j,i) = P(n,k,j,ie-ngh+1+i);
	  Bs.Xe(n,k,j,i) = P(n,k,j,is      +i);
	  //Bs.Xs_data[((n*ktot+k)*jtot+j)*ngh+i]=P.data[((n*ktot+k)*jtot+j)*itot+ie-ngh+1+i];
	  //Bs.Xe_data[((n*ktot+k)*jtot+j)*ngh+i]=P.data[((n*ktot+k)*jtot+j)*itot+is      +i];
  }
  
  // y-direction
  // |     |Bs.Ye   Bs.Ys|     |
  // |Br.Ys|             |Br.Ye|
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=0; j<ngh;j++)
	for (int i=is; i<=ie; i++) {
	  Bs.Ys(n,k,j,i) = P(n,k,je-ngh+1+j,i);
	  Bs.Ye(n,k,j,i) = P(n,k,js      +j,i);
	  //Bs.Ys_data[((n*ktot+k)*ngh+j)*itot+i]=P.data[((n*ktot+k)*jtot+je-ngh+1+j)*itot+i];
	  //Bs.Ye_data[((n*ktot+k)*ngh+j)*itot+i]=P.data[((n*ktot+k)*jtot+js      +j)*itot+i];
  }
  
  // z-direction
  // |     |Bs.Ze   Bs.Zs|     |
  // |Br.Zs|             |Br.Ze|
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  Bs.Zs(n,k,j,i) = P(n,ke-ngh+1+k,j,i);
	  Bs.Ze(n,k,j,i) = P(n,ks      +k,j,i);
	  //Bs.Zs_data[((n*ngh+k)*jtot+j)*itot+i]=P.data[((n*ktot+ke-ngh+1+k)*jtot+j)*itot+i];
	  //Bs.Ze_data[((n*ngh+k)*jtot+j)*itot+i]=P.data[((n*ktot+ks      +k)*jtot+j)*itot+i];
  }

  SendRecvBoundary(Bs,Br);

  
  // |     |Bs.Xe   Bs.Xs|     |
  // |Br.Xs|             |Br.Xe|
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

  // |     |Bs.Ye   Bs.Ys|     |
  // |Br.Ys|             |Br.Ye|
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

  // |     |Bs.Ze   Bs.Zs|     |
  // |Br.Zs|             |Br.Ze|
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

};// end namespace
