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
  Array4D<double> Fx,Fy,Fz;
  int nprim{5}; //!
  Array4D<double> P; //! P(nprim,ktot,jtot,itot)
  int nden{0},nvex{1},nvey{2},nvez{3},nene{4};
#pragma acc declare create ( U)
#pragma acc declare create ( P)
#pragma acc declare create (Fx)
#pragma acc declare create (Fy)
#pragma acc declare create (Fz)

};

using namespace hydflux_mod;

void AllocateVariables(){
  U.allocate(mconsv,ktot,jtot,itot);
  P.allocate(nprim ,ktot,jtot,itot);
      
  Fx.allocate(mconsv,ktot,jtot,itot);
  Fy.allocate(mconsv,ktot,jtot,itot);
  Fz.allocate(mconsv,ktot,jtot,itot);


  for (int m=0; m<mconsv; m++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  U(m,k,j,i) = 0.0;
	  Fx(m,k,j,i) = 0.0;
	  Fy(m,k,j,i) = 0.0;
	  Fz(m,k,j,i) = 0.0;
  }

  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  P(n,k,j,i) = 0.0;
  }
  
#pragma acc update device( U.data[0: U.size], U.n1, U.n2, U.n3, U.nv)
#pragma acc update device(Fx.data[0:Fx.size],Fx.n1,Fx.n2,Fx.n3,Fx.nv)
#pragma acc update device(Fy.data[0:Fy.size],Fy.n1,Fy.n2,Fy.n3,Fy.nv)
#pragma acc update device(Fz.data[0:Fz.size],Fz.n1,Fz.n2,Fz.n3,Fz.nv)
#pragma acc update device( P.data[0: P.size], P.n1, P.n2, P.n3, P.nv)
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
	Fx(mden,k,j,i) = ux * q_up;
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
	Fy(mden,k,j,i) = uy * q_up;
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
	Fz(mden,k,j,i) = uz * q_up;
      }
}

void UpdateConservU(){

#pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
        double dqx = (Fx(mden,k,j,i+1) - Fx(mden,k,j,i)) / dx;
        double dqy = (Fy(mden,k,j+1,i) - Fy(mden,k,j,i)) / dy;
        double dqz = (Fz(mden,k+1,j,i) - Fz(mden,k,j,i)) / dz;
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
