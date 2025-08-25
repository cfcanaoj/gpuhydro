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

  inline constexpr int mconsv{9},madd{3}; //!
  inline constexpr int mudn{ 0},muvu{ 1},muvv{ 2},muvw{ 3},muet{ 4},
                                mubu{ 5},mubv{ 6},mubw{ 7},mubp{ 8},
                       mfdn{ 9},mfvu{10},mfvv{11},mfvw{12},mfet{13},
                       mfbu{14},mfbv{15},mfbw{16},mfbp{17},
                       mcsp{18},mvel{19},mpre{20};
  inline constexpr int mden{ 0},mrv1{ 1},mrv2{ 2},mrv3{ 3},meto{ 4},
                                mbm1{ 5},mbm2{ 6},mbm3{ 7},mbps{ 8},
                                mrvu{ 1},mrvv{ 2},mrvw{ 3},
                                mbmu{ 5},mbmv{ 6},mbmw{ 7};
  inline constexpr int nprim{11}; //!
  inline constexpr int nden{0},nve1{1},nve2{2},nve3{3},nene{4},npre{5},ncsp{6},
                               nbm1{7},nbm2{8},nbm3{9},nbps{10};

  extern Array4D<double> P; //! P(nprim ,ktot,jtot,itot)
  extern Array4D<double> U; //! U(mconsv,ktot,jtot,itot)
  extern Array4D<double> Fx,Fy,Fz;
  extern double csiso;
  extern double chg;
  
};

void AllocateHydroVariables(Array4D<double>& U, Array4D<double>& Fx,Array4D<double>& Fy,Array4D<double>& Fz,Array4D<double>& P);
void DeallocateHydroVariables(Array4D<double>& U, Array4D<double>& Fx,Array4D<double>& Fy,Array4D<double>& Fz,Array4D<double>& P);
void GetNumericalFlux1(const Array4D<double>& P,Array4D<double>& Fx);
void GetNumericalFlux2(const Array4D<double>& P,Array4D<double>& Fy);
void GetNumericalFlux3(const Array4D<double>& P,Array4D<double>& Fz);
void UpdateConservU(const Array4D<double>& Fx,const Array4D<double>& Fy,const Array4D<double>& Fz,Array4D<double>& U);
void UpdatePrimitvP(const Array4D<double>& U,Array4D<double>& P);
void ControlTimestep();
void EvaluateCh();
void DampPsi(Array4D<double>& U);

#endif
