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
  extern BoundaryArray<double> Bs,Br; 
};

void AllocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br);
void DeallocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br);

void SetBoundaryCondition(FieldArray<double>& P,BoundaryArray<double>& Bs,BoundaryArray<double>& Br);

void SendRecvBoundary(const BoundaryArray<double>& Bs,BoundaryArray<double>& Br);

#endif
