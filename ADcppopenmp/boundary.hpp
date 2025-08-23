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
  extern Array4D<double> Xs,Xe,Ys,Ye,Zs,Ze; 
};

void AllocateBoundaryVariables(Array4D<double>& Xs,Array4D<double>& Xe
			      ,Array4D<double>& Ys,Array4D<double>& Ye
			      ,Array4D<double>& Zs,Array4D<double>& Ze);
void DeallocateBoundaryVariables(Array4D<double>& Xs,Array4D<double>& Xe
			      ,Array4D<double>& Ys,Array4D<double>& Ye
			      ,Array4D<double>& Zs,Array4D<double>& Ze);

void SetBoundaryCondition(Array4D<double>& P
			 ,Array4D<double>& Xs,Array4D<double>& Xe
			 ,Array4D<double>& Ys,Array4D<double>& Ye
			 ,Array4D<double>& Zs,Array4D<double>& Ze);

#endif
