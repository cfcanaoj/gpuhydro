/**
 * @file boundary.hpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/
#ifndef BOUNDARY_HPP_
#define BOUNDARY_HPP_

#include "hydro.hpp"
using namespace hydflux_mod;

namespace boundary_mod {

  template <typename T>
  class BoundaryArray {
  public:
    T* Xs_data = nullptr;
    T* Xe_data = nullptr;
    T* Ys_data = nullptr;
    T* Ye_data = nullptr;
    T* Zs_data = nullptr;
    T* Ze_data = nullptr;
    int nv=0, ng=0, n3 = 0, n2 = 0, n1 = 0;
    size_t size1 = 0;
    size_t size2 = 0;
    size_t size3 = 0;
    
    BoundaryArray() = default;
    BoundaryArray(int nv_,int ng_,int n3_, int n2_, int n1_) {
      allocate(nv_,ng_,n3_,n2_,n1_);
    }
    void allocate(int nv_,int ng_,int n3_, int n2_, int n1_) {
      nv = nv_;ng = ng_; n3 = n3_; n2 = n2_; n1 = n1_;
      size1 = static_cast<size_t>(nv) * n3 * n2 * ng;
      Xs_data = static_cast<T*>(malloc(sizeof(T) * size1));
      Xe_data = static_cast<T*>(malloc(sizeof(T) * size1));
      size2 = static_cast<size_t>(nv) * n3 * ng * n1;
      Ys_data = static_cast<T*>(malloc(sizeof(T) * size2));
      Ye_data = static_cast<T*>(malloc(sizeof(T) * size2));
      size3 = static_cast<size_t>(nv) * ng * n2 * n1;
      Zs_data = static_cast<T*>(malloc(sizeof(T) * size3));
      Ze_data = static_cast<T*>(malloc(sizeof(T) * size3));
    }
    
    inline       T& Xs(int n, int k, int j, int i)       noexcept { return Xs_data[((n*n3 + k)*n2 + j)*ng + i]; }
    inline const T& Xs(int n, int k, int j, int i) const noexcept { return Xs_data[((n*n3 + k)*n2 + j)*ng + i]; }
    inline       T& Xe(int n, int k, int j, int i)       noexcept { return Xe_data[((n*n3 + k)*n2 + j)*ng + i]; }
    inline const T& Xe(int n, int k, int j, int i) const noexcept { return Xe_data[((n*n3 + k)*n2 + j)*ng + i]; }
    inline       T& Ys(int n, int k, int j, int i)       noexcept { return Ys_data[((n*n3 + k)*ng + j)*n1 + i]; }
    inline const T& Ys(int n, int k, int j, int i) const noexcept { return Ys_data[((n*n3 + k)*ng + j)*n1 + i]; }
    inline       T& Ye(int n, int k, int j, int i)       noexcept { return Ye_data[((n*n3 + k)*ng + j)*n1 + i]; }
    inline const T& Ye(int n, int k, int j, int i) const noexcept { return Ye_data[((n*n3 + k)*ng + j)*n1 + i]; }
    inline       T& Zs(int n, int k, int j, int i)       noexcept { return Zs_data[((n*ng + k)*n2 + j)*n1 + i]; }
    inline const T& Zs(int n, int k, int j, int i) const noexcept { return Zs_data[((n*ng + k)*n2 + j)*n1 + i]; }
    inline       T& Ze(int n, int k, int j, int i)       noexcept { return Ze_data[((n*ng + k)*n2 + j)*n1 + i]; }
    inline const T& Ze(int n, int k, int j, int i) const noexcept { return Ze_data[((n*ng + k)*n2 + j)*n1 + i]; }
    
  };

  extern BoundaryArray<double> Bs,Br; 

  void AllocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br);
  void DeallocateBoundaryVariables(BoundaryArray<double>& Bs,BoundaryArray<double>& Br);
  
  void SetBoundaryCondition(FieldArray<double>& P,BoundaryArray<double>& Bs,BoundaryArray<double>& Br);

  void SendRecvBoundary(const BoundaryArray<double>& Bs,BoundaryArray<double>& Br);
};

#endif
