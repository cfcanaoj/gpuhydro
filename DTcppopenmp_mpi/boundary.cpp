
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
#include "hydro_arrays.hpp"
#include "resolution.hpp"
#include "hydro.hpp"

using namespace hydro_arrays_mod;

namespace boundary_mod {
#pragma omp declare target
  Array4D<double> Xs,Xe,Ys,Ye,Zs,Ze;
#pragma omp end declare target 
};

using namespace hydflux_mod;

void AllocateBoundaryVariables(Array4D<double>& Xs,Array4D<double>& Xe
			      ,Array4D<double>& Ys,Array4D<double>& Ye
			      ,Array4D<double>& Zs,Array4D<double>& Ze){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  Xs.allocate(nprim ,ktot,jtot,ngh );
  Xe.allocate(nprim ,ktot,jtot,ngh );
  Ys.allocate(nprim ,ktot,ngh ,itot);
  Ye.allocate(nprim ,ktot,ngh ,itot);
  Zs.allocate(nprim ,ngh ,jtot,itot);
  Ze.allocate(nprim ,ngh ,jtot,itot);


#pragma omp target enter data map (alloc: Xs.data[0: Xs.size])
#pragma omp target enter data map (alloc: Xe.data[0: Xe.size])
#pragma omp target enter data map (alloc: Ys.data[0: Ys.size])
#pragma omp target enter data map (alloc: Ye.data[0: Ye.size])
#pragma omp target enter data map (alloc: Zs.data[0: Zs.size])
#pragma omp target enter data map (alloc: Ze.data[0: Ze.size])  
#pragma omp target update to (Xs.data[0: Xs.size], Xs.n1, Xs.n2, Xs.n3, Xs.nv)
#pragma omp target update to (Xe.data[0: Xe.size], Xe.n1, Xe.n2, Xe.n3, Xe.nv)
#pragma omp target update to (Ys.data[0: Ys.size], Ys.n1, Ys.n2, Ys.n3, Ys.nv)
#pragma omp target update to (Ye.data[0: Ye.size], Ye.n1, Ye.n2, Ye.n3, Ye.nv)
#pragma omp target update to (Zs.data[0: Zs.size], Zs.n1, Zs.n2, Zs.n3, Zs.nv)
#pragma omp target update to (Ze.data[0: Ze.size], Ze.n1, Ze.n2, Ze.n3, Ze.nv)  

}


void DeallocateBoundaryVariables(Array4D<double>& Xs,Array4D<double>& Xe
			      ,Array4D<double>& Ys,Array4D<double>& Ye
			      ,Array4D<double>& Zs,Array4D<double>& Ze){
  using namespace resolution_mod;
  using namespace hydflux_mod;


#pragma omp target exit data map (delete: Xs.data[0: Xs.size], Xs.n1, Xs.n2, Xs.n3, Xs.nv)
#pragma omp target exit data map (delete: Xe.data[0: Xe.size], Xe.n1, Xe.n2, Xe.n3, Xe.nv)
#pragma omp target exit data map (delete: Ys.data[0: Ys.size], Ys.n1, Ys.n2, Ys.n3, Ys.nv)
#pragma omp target exit data map (delete: Ye.data[0: Ye.size], Ye.n1, Ye.n2, Ye.n3, Ye.nv)
#pragma omp target exit data map (delete: Zs.data[0: Zs.size], Zs.n1, Zs.n2, Zs.n3, Zs.nv)
#pragma omp target exit data map (delete: Ze.data[0: Ze.size], Ze.n1, Ze.n2, Ze.n3, Ze.nv)  
}



void SetBoundaryCondition(Array4D<double>& P,Array4D<double>& Xs,Array4D<double>& Xe
			                    ,Array4D<double>& Ys,Array4D<double>& Ye
			                    ,Array4D<double>& Zs,Array4D<double>& Ze) {
  using namespace resolution_mod;
  using namespace hydflux_mod;

#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=0; i<ngh; i++) {
	  //	  Xs(n,k,j,i) = P(n,k,j,ie-ngh+1+i);
	  //Xe(n,k,j,i) = P(n,k,j,is      +i);
	  Xs.data[((n*ktot+k)*jtot+j)*ngh+i]=P.data[((n*ktot+k)*jtot+j)*itot+ie-ngh+1+i];
	  Xe.data[((n*ktot+k)*jtot+j)*ngh+i]=P.data[((n*ktot+k)*jtot+j)*itot+is      +i];
  }
  
  // y-direction
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=0; j<ngh;j++)
	for (int i=is; i<=ie; i++) {
	  //Ys(n,k,j,i) = P(n,k,je-ngh+1+j,i);
	  //Ye(n,k,j,i) = P(n,k,js      +j,i);
	  Ys.data[((n*ktot+k)*ngh+j)*itot+i]=P.data[((n*ktot+k)*jtot+je-ngh+1+j)*itot+i];
	  Ye.data[((n*ktot+k)*ngh+j)*itot+i]=P.data[((n*ktot+k)*jtot+js      +j)*itot+i];
  }
  
  // z-direction
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  //Zs(n,k,j,i) = P(n,ke-ngh+1+k,j,i);
	  //Ze(n,k,j,i) = P(n,ks      +k,j,i);
	  Zs.data[((n*ngh+k)*jtot+j)*itot+i]=P.data[((n*ktot+ke-ngh+1+k)*jtot+j)*itot+i];
	  Ze.data[((n*ngh+k)*jtot+j)*itot+i]=P.data[((n*ktot+ks      +k)*jtot+j)*itot+i];
  }


#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=0; i<ngh; i++) {
	  //P(n,k,j,is-ngh+i) = Xs(n,k,j,i);
	  //P(n,k,j,ie+1  +i) = Xe(n,k,j,i);
	  P.data[((n*ktot+k)*jtot+j)*itot+is-ngh+i] = Xs.data[((n*ktot+k)*jtot+j)*ngh+i];
	  P.data[((n*ktot+k)*jtot+j)*itot+ie+1  +i] = Xe.data[((n*ktot+k)*jtot+j)*ngh+i];
  }

#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=0; j<ngh;j++)
	for (int i=is; i<=ie; i++) {
	  //P(n,k,js-ngh+j,i) = Ys(n,k,j,i);
	  //P(n,k,je+1  +j,i) = Ye(n,k,j,i);
	  P.data[((n*ktot+k)*jtot+js-ngh+j)*itot+i] = Ys.data[((n*ktot+k)*ngh+j)*itot+i];
	  P.data[((n*ktot+k)*jtot+je+1  +j)*itot+i] = Ye.data[((n*ktot+k)*ngh+j)*itot+i];
  }

#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  //P(n,ks-ngh+k,j,i) = Zs(n,k,j,i);
	  //P(n,ke+1  +k,j,i) = Ze(n,k,j,i);
	  P.data[((n*ktot+ks-ngh+k)*jtot+j)*itot+i] = Zs.data[((n*ngh+k)*jtot+j)*itot+i];
	  P.data[((n*ktot+ke+1  +k)*jtot+j)*itot+i] = Ze.data[((n*ngh+k)*jtot+j)*itot+i];
  }
};

