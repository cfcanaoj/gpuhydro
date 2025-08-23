/**
 * @file hydro.cpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/

#include <cmath>
#include <cstdio>
#include <algorithm>
#include "hydro_arrays.hpp"
#include "resolution.hpp"
#include "hydro.hpp"

using namespace hydro_arrays_mod;
using namespace resolution_mod;

namespace hydflux_mod {
  int mconsv{5}; //!
  Array4D<double> U; //! U(mconsv,ktot,jtot,itot)
  int mden{0},mrvx{1},mrvy{2},mrvz{3},meto{4};
  Array4D<double> fluxx,fluxy,fluxz;
  int nprim{5}; //!
  Array4D<double> P; //! P(nprim,ktot,jtot,itot)
  int nden{0},nvex{1},nvey{2},nvez{3},nene{4};
#pragma acc declare create (U,P)
#pragma acc declare create (fluxx,fluxy,fluxz)

};

using namespace hydflux_mod;

void AllocateVariables(){
  U.allocate(mconsv,ktot,jtot,itot);
  P.allocate(nprim ,ktot,jtot,itot);
      
  fluxx.allocate(mconsv,ktot,jtot,itot);
  fluxy.allocate(mconsv,ktot,jtot,itot);
  fluxz.allocate(mconsv,ktot,jtot,itot);

  printf("alloc1\n");
  for (int m=0; m<mconsv; m++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  U(m,k,j,i) = 0.0;
	  fluxx(m,k,j,i) = 0.0;
	  fluxy(m,k,j,i) = 0.0;
	  fluxz(m,k,j,i) = 0.0;
  }
  printf("alloc2\n");
  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  P(n,k,j,i) = 0.0;
  }
  
#pragma acc update device(U.data[0:fluxx.size()],U.n1,U.n2,U.n3,U.nv)
#pragma acc update device(fluxx.data[0:fluxx.size()],fluxx.n1,fluxx.n2,fluxx.n3,fluxx.nv)
#pragma acc update device(fluxy.data[0:fluxy.size()],fluxy.n1,fluxy.n2,fluxy.n3,fluxy.nv)
#pragma acc update device(fluxz.data[0:fluxz.size()],fluxz.n1,fluxz.n2,fluxz.n3,fluxz.nv)
#pragma acc update device(P.data[0:fluxx.size()],P.n1,P.n2,P.n3,P.nv)
}

void GetNumericalFlux1(){

#pragma acc loop independent  collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie+1; ++i) {
	double qL = P(nden,k,j,i-1);
	double qR = P(nden,k,j,i  );
	double ux = 0.5e0*(P(nvex,k,j,i-1)+P(nvex,k,j,i));
	double q_up = (ux >= 0.0) ? qL : qR;
	fluxx(mden,k,j,i) = ux * q_up;
      }
}


void GetNumericalFlux2(){

#pragma acc loop independent collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je+1; ++j)
      for (int i=is; i<=ie; ++i) {
	double qL = P(nden,k,j-1,i);
	double qR = P(nden,k,j  ,i);
	double uy = 0.5e0*(P(nvey,k,j-1,i)+P(nvey,k,j,i));
	double q_up = (uy >= 0.0) ? qL : qR;
	fluxy(mden,k,j,i) = uy * q_up;
      }
}

void GetNumericalFlux3(){

#pragma acc loop independent collapse(3)
  for (int k=ks; k<=ke+1; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
	double qL = P(nden,k-1,j,i);
	double qR = P(nden,k  ,j,i);
	double uz = 0.5e0*(P(nvez,k-1,j,i)+P(nvez,k,j,i));
	double q_up = (uz >= 0.0) ? qL : qR;
	fluxz(mden,k,j,i) = uz * q_up;
      }
}

void UpdateConservU(){

#pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
        double dqx = (fluxx(mden,k,j,i+1) - fluxx(mden,k,j,i)) / dx;
        double dqy = (fluxy(mden,k,j+1,i) - fluxy(mden,k,j,i)) / dy;
        double dqz = (fluxz(mden,k+1,j,i) - fluxz(mden,k,j,i)) / dz;
        U(mden,k,j,i) -= dt * (dqx + dqy + dqz);
      }
}


void UpdatePrimitvP(){
#pragma acc loop independent collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
         P(nden,k,j,i) = U(mden,k,j,i);
	 //P(nvex,k,j,i) = U(mrvx,k,j,i)/U(mden,k,j,i);
         //P(nvey,k,j,i) = U(mrvy,k,j,i)/U(mden,k,j,i);
         //P(nvez,k,j,i) = U(mrvz,k,j,i)/U(mden,k,j,i);
	 double ekin = 0.5e0*( U(mrvx,k,j,i)*U(mrvx,k,j,i)
			      +U(mrvy,k,j,i)*U(mrvy,k,j,i)
			      +U(mrvz,k,j,i)*U(mrvz,k,j,i))/U(mden,k,j,i);
         P(nene,k,j,i) = U(meto,k,j,i)-ekin;
      }
}

void ControlTimestep(){
  double dtmin;
  const double eps = 1.0e-10;
  dtmin = 1.0e10;
  for (int k=ks; k<=ke; k++)
    for (int j=js; j<=je; j++)
      for (int i=is; i<=ie; i++) {
	double dtminloc = std::min({dx/(std::abs(P(nvex,k,j,i))+eps)
				    ,dy/(std::abs(P(nvey,k,j,i))+eps)
				   ,dz/(std::abs(P(nvez,k,j,i))+eps)});
	dtmin = std::min(dtminloc,dtmin);
      }
  dt = 0.5*dtmin;
}
