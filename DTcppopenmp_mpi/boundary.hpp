/**
 * @file boundary.hpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/
#ifndef BOUNDARY_HPP_
#define BOUNDARY_HPP_

#include "hydro_arrays.hpp"
using namespace hydro_arrays_mod;

namespace boundary_mod {
  extern Boundary3D<double> Bs,Br; 
};

void AllocateBoundaryVariables(Boundary3D<double>& Bs,Boundary3D<double>& Br);
void DeallocateBoundaryVariables(Boundary3D<double>& Bs,Boundary3D<double>& Br);

void SetBoundaryCondition(Array4D<double>& P,Boundary3D<double>& Bs,Boundary3D<double>& Br);

void SendRecvBoundary(const Boundary3D<double>& Bs,Boundary3D<double>& Br);

#endif
