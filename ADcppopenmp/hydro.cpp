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
#pragma omp declare target
  int mconsv{5}; //!
  HydroArrays<double> U; //! U(mconsv,ktot,jtot,itot)
  int mden{0},mrvx{1},mrvy{2},mrvz{3},meto{4};
  HydroArrays<double> fluxx,fluxy,fluxz;
  int nprim{5}; //!
  HydroArrays<double> P; //! P(nprim,ktot,jtot,itot)
  int nden{0},nvex{1},nvey{2},nvez{3},nene{4};
#pragma omp end declare target 
};

using namespace hydflux_mod;

void AllocateVariables(){
  U.allocate(mconsv,ktot,jtot,itot);
  P.allocate(nprim ,ktot,jtot,itot);
      
  fluxx.allocate(mconsv,ktot,jtot,itot);
  fluxy.allocate(mconsv,ktot,jtot,itot);
  fluxz.allocate(mconsv,ktot,jtot,itot);

#pragma omp target update to (fluxx.data()[0:fluxx.size()])
#pragma omp target update to (fluxy.data()[0:fluxy.size()])
#pragma omp target update to (fluxz.data()[0:fluxz.size()])

}

void GetNumericalFlux1(){

#pragma omp target teams distribute parallel for collapse(3)
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

#pragma omp target teams distribute parallel for collapse(3)
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

#pragma omp target teams distribute parallel for collapse(3)
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

#pragma omp target teams distribute parallel for collapse(3)
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
