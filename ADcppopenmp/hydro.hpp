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

namespace hydflux_mod {
  using namespace hydro_arrays_mod;
  extern int mconsv; //!
  extern Array4D<double> U; //! U(mconsv,ktot,jtot,itot)
  extern int mden,mrvx,mrvy,mrvz,meto;
  extern Array4D<double> Fx,Fy,Fz;
  extern int nprim; //!
  extern Array4D<double> P; //! P(nprim,ktot,jtot,itot)
  extern int nden,nvex,nvey,nvez,nene;
};

void AllocateVariables();
void GetNumericalFlux1();
void GetNumericalFlux2();
void GetNumericalFlux3();
void UpdateConservU();
void UpdatePrimitvP();
void ControlTimestep();

#endif
