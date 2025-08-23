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
  Array4D<double> Xs,Xe,Ys,Ye,Zs,Ze; 
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

#pragma omp target update to ( Xs.data[0: Xs.size()], Xs.n1, Xs.n2, Xs.n3, Xs.nv)
#pragma omp target update to ( Xe.data[0: Xe.size()], Xe.n1, Xe.n2, Xe.n3, Xe.nv)
#pragma omp target update to ( Ys.data[0: Ys.size()], Ys.n1, Ys.n2, Ys.n3, Ys.nv)
#pragma omp target update to ( Ye.data[0: Ye.size()], Ye.n1, Ye.n2, Ye.n3, Ye.nv)
#pragma omp target update to ( Zs.data[0: Zs.size()], Zs.n1, Zs.n2, Zs.n3, Zs.nv)
#pragma omp target update to ( Ze.data[0: Ze.size()], Ze.n1, Ze.n2, Ze.n3, Ze.nv)  
}

void SetBoundaryCondition(Array4D<double>& P,Array4D<double>& Xs,Array4D<double>& Xe
			                    ,Array4D<double>& Ys,Array4D<double>& Ye
			                    ,Array4D<double>& Zs,Array4D<double>& Ze) {
  using namespace resolution_mod;
  using namespace hydflux_mod;
  
  // x-direction
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=0; i<ngh; i++) {
	  Xs(n,k,j,i) = P(n,k,j,ie-ngh+i);
	  Xe(n,k,j,i) = P(n,k,j,is    +i);
  }
  // y-direction
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=0; j<ngh;j++)
	for (int i=is; i<=ie; i++) {
	  Ys(n,k,j,i) = P(n,k,je-ngh+j,i);
	  Ye(n,k,j,i) = P(n,k,js    +j,i);
  }
  
  // z-direction
#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  Zs(n,k,j,i) = P(n,ke-ngh+k,j,i);
	  Ze(n,k,j,i) = P(n,ks    +k,j,i);
  }


#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=0; i<ngh; i++) {
	  P(n,k,j,is-ngh+i) = Xs(n,k,j,i);
	  P(n,k,j,ie+1  +i) = Xe(n,k,j,i);
  }

#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=0; j<ngh;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,k,js-ngh+j,i) = Ys(n,k,j,i);
	  P(n,k,je+1  +j,i) = Ye(n,k,j,i);
  }

#pragma omp target teams distribute parallel for collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,ks-ngh+k,j,i) = Zs(n,k,j,i);
	  P(n,ke+1  +k,j,i) = Ze(n,k,j,i);
  }


}

