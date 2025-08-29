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

  extern GridArray<double> G;
  extern FieldArray<double> P; //! P(nprim ,ktot,jtot,itot)
  extern FieldArray<double> U; //! U(mconsv,ktot,jtot,itot)
  extern FieldArray<double> Fx,Fy,Fz;
  extern double csiso;
  extern double chg;
  
};

void AllocateHydroVariables(GridArray<double>& G,FieldArray<double>& U,FieldArray<double>& Fx,FieldArray<double>& Fy,FieldArray<double>& Fz,FieldArray<double>& P);
void DeallocateHydroVariables(GridArray<double> G,FieldArray<double>& U, FieldArray<double>& Fx,FieldArray<double>& Fy,FieldArray<double>& Fz,FieldArray<double>& P);
void GetNumericalFlux1(const GridArray<double>& G,const FieldArray<double>& P,FieldArray<double>& Fx);
void GetNumericalFlux2(const GridArray<double>& G,const FieldArray<double>& P,FieldArray<double>& Fy);
void GetNumericalFlux3(const GridArray<double>& G,const FieldArray<double>& P,FieldArray<double>& Fz);
void UpdateConservU(const GridArray<double>& G,const FieldArray<double>& Fx,const FieldArray<double>& Fy,const FieldArray<double>& Fz,FieldArray<double>& U);
void UpdatePrimitvP(const FieldArray<double>& U,FieldArray<double>& P);
void ControlTimestep(const GridArray<double>& G);
void EvaluateCh();
void DampPsi(const GridArray<double>& G,FieldArray<double>& U);

#endif
