/**
 * @file hydro.cpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/
#ifndef HYDRO_HPP_
#define HYDRO_HPP_

#include <cmath>
#include <cstdio>
#include <algorithm>
#include "hydro_arrays.hpp"

using namespace hydro_arrays_mod;

namespace hydflux_mod {
  extern int mconsv; //!
  extern Array4D<double> U; //! U(mconsv,ktot,jtot,itot)
  extern int mden,mrvx,mrvy,mrvz,meto;
  extern Array4D<double> Fx,Fy,Fz;
  extern int nprim; //!
  extern Array4D<double> P; //! P(nprim,ktot,jtot,itot)
  extern int nden,nvex,nvey,nvez,nene;
};

void AllocateHydroVariables(Array4D<double>& U, Array4D<double>& Fx,Array4D<double>& Fy,Array4D<double>& Fz,Array4D<double>& P);
void DeallocateHydroVariables(Array4D<double>& U, Array4D<double>& Fx,Array4D<double>& Fy,Array4D<double>& Fz,Array4D<double>& P);
void GetNumericalFlux1(const Array4D<double>& P,Array4D<double>& Fx);
void GetNumericalFlux2(const Array4D<double>& P,Array4D<double>& Fy);
void GetNumericalFlux3(const Array4D<double>& P,Array4D<double>& Fz);
void UpdateConservU(const Array4D<double>& Fx,const Array4D<double>& Fy,const Array4D<double>& Fz,Array4D<double>& U);
void UpdatePrimitvP(const Array4D<double>& U,Array4D<double>& P);
void ControlTimestep();

#endif
