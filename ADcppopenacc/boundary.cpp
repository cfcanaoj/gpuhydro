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

void SetBoundaryCondition() {
  using namespace resolution_mod;
  using namespace hydflux_mod;
  double dummy;
  
  // x-direction
#pragma acc loop independent collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=1; i<=ngh; i++) {
	  P(n,k,j,is-i) = P(n,k,j,ie+1-i);
	  P(n,k,j,ie+i) = P(n,k,j,is-1+i);
	  //dummy = P(n,k,j,ie+1-i);
	  //dummy = 1.0;
  }
  // y-direction
#pragma acc loop independent collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=1; j<=ngh;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,k,js-j,i) = P(n,k,je+1-j,i);
	  P(n,k,je+j,i) = P(n,k,js-1+j,i);
  }
  
  // z-direction
#pragma acc loop independent collapse(4)
  for (int n=0; n<nprim; n++)
    for (int k=1; k<=ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,ks-k,j,i) = P(n,ke+1-k,j,i);
	  P(n,ke+k,j,i) = P(n,ks-1+k,j,i);
  }
  
}

