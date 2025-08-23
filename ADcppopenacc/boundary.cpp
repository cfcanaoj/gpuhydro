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
using namespace resolution_mod;

namespace boundary_mod{

  Array4D<double> Xs,Xe;
  Array4D<double> Ys,Ye;
  Array4D<double> Zs,Ze;
#pragma acc declare create (Xs,Xe)
#pragma acc declare create (Ys,Ye)
#pragma acc declare create (Zs,Ze)

}

using namespace boundary_mod;

void SetBoundaryCondition() {
  using namespace resolution_mod;
  using namespace hydflux_mod;
  static bool is_inited= false;

  if(! is_inited){   
    Xs.allocate(nprim ,ktot,jtot,ngh );
    Xe.allocate(nprim ,ktot,jtot,ngh );
    Ys.allocate(nprim ,ktot,ngh ,itot);
    Ye.allocate(nprim ,ktot,ngh ,itot);
    Zs.allocate(nprim ,ngh ,jtot,itot);
    Ze.allocate(nprim ,ngh ,jtot,itot);
    is_inited = true;
#pragma acc update device(Xs.data[0:Xs.size()],Xs.n1,Xs.n2,Xs.n3,Xs.nv)
#pragma acc update device(Xe.data[0:Xe.size()],Xe.n1,Xe.n2,Xe.n3,Xe.nv)
#pragma acc update device(Ys.data[0:Ys.size()],Ys.n1,Ys.n2,Ys.n3,Ys.nv)
#pragma acc update device(Ye.data[0:Ye.size()],Ye.n1,Ye.n2,Ye.n3,Ye.nv)
#pragma acc update device(Zs.data[0:Zs.size()],Zs.n1,Zs.n2,Zs.n3,Zs.nv)
#pragma acc update device(Ze.data[0:Ze.size()],Ze.n1,Ze.n2,Ze.n3,Ze.nv)
  }
  
#pragma acc data present(P.data[0:P.size()])
  {
  // x-direction
#pragma acc loop independent collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=0; i<ngh; i++) {
	  Xs(n,k,j,i) = P(n,k,j,ie-ngh+i);
	  Xe(n,k,j,i) = P(n,k,j,is    +i);
  }
  // y-direction
#pragma acc loop independent collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=0; j<ngh;j++)
	for (int i=is; i<=ie; i++) {
	  Ys(n,k,j,i) = P(n,k,je-ngh+j,i);
	  Ye(n,k,j,i) = P(n,k,js    +j,i);
  }
  
  // z-direction
#pragma acc loop independent collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  Zs(n,k,j,i) = P(n,ke-ngh+k,j,i);
	  Ze(n,k,j,i) = P(n,ks    +k,j,i);
  }


#pragma acc loop independent collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=0; i<ngh; i++) {
	  P(n,k,j,is-ngh+i) = Xs(n,k,j,i);
	  P(n,k,j,ie+1  +i) = Xe(n,k,j,i);
  }

#pragma acc loop independent collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=0; j<ngh;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,k,js-ngh+j,i) = Ys(n,k,j,i);
	  P(n,k,je+1  +j,i) = Ye(n,k,j,i);
  }

#pragma acc loop independent collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,ks-ngh+k,j,i) = Zs(n,k,j,i);
	  P(n,ke+1  +k,j,i) = Ze(n,k,j,i);
  }


  }//! acc end data
}

