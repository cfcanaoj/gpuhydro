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

namespace hydflux_mod {
  using index_t = int;     // 必要なら 32/64 を切り替え
  using size_t  = std::size_t;
  
  template <typename T>
  class GridArray {
  public:
    T* x1a_data = nullptr;
    T* x1b_data = nullptr;
    T* x2a_data = nullptr;
    T* x2b_data = nullptr;
    T* x3a_data = nullptr;
    T* x3b_data = nullptr;
    int n3 = 0, n2 = 0, n1 = 0;
    size_t size = 0;
    GridArray() = default;
    GridArray(int n3_, int n2_, int n1_) {
      allocate(n3_,n2_,n1_);
    }
    void allocate(int n3_, int n2_, int n1_) {
      n3 = n3_; n2 = n2_; n1 = n1_;
      size = static_cast<size_t>(n1);
      x1a_data = static_cast<T*>(malloc(sizeof(T) * size));
      x1b_data = static_cast<T*>(malloc(sizeof(T) * size));
      size = static_cast<size_t>(n2);
      x2a_data = static_cast<T*>(malloc(sizeof(T) * size));
      x2b_data = static_cast<T*>(malloc(sizeof(T) * size));
      size = static_cast<size_t>(n3);
      x3a_data = static_cast<T*>(malloc(sizeof(T) * size));
      x3b_data = static_cast<T*>(malloc(sizeof(T) * size));
    }
    
    inline T& x1a(int i)       noexcept { return x1a_data[i]; }
    inline T& x1b(int i)       noexcept { return x1b_data[i]; }
    // const
    inline const T& x1a(int i) const noexcept { return x1a_data[i]; }
    inline const T& x1b(int i) const noexcept { return x1b_data[i]; }

    inline T& x2a(int j)       noexcept { return x2a_data[j]; }
    inline T& x2b(int j)       noexcept { return x2b_data[j]; }
    // const
    inline const T& x2a(int j) const noexcept { return x2a_data[j]; }
    inline const T& x2b(int j) const noexcept { return x2b_data[j]; }
  
    inline T& x3a(int k)       noexcept { return x3a_data[k]; }
    inline T& x3b(int k)       noexcept { return x3b_data[k]; }
    // const
    inline const T& x3a(int k) const noexcept { return x3a_data[k]; }
    inline const T& x3b(int k) const noexcept { return x3b_data[k]; }

  };

  
  template <typename T>
  class FieldArray {
  public:
    T* data = nullptr;
    int nv = 0, n3 = 0, n2 = 0, n1 = 0;
    size_t size = 0;
    FieldArray() = default;
    FieldArray(int _nv, int _n3, int _n2, int _n1) {
      allocate(_nv,_n3,_n2,_n1);
    }
    
    void allocate(int _nv, int _n3, int _n2, int _n1) {
      nv = _nv; n3 = _n3; n2 = _n2; n1 = _n1;
      size = static_cast<size_t>(nv) * n3 * n2 * n1;
      data = static_cast<T*>(malloc(sizeof(T) * size));
  }
    
    void deallocate() {
      free(data);
      data = nullptr;
    }
    
    inline const T& operator()(int n, int k, int j, int i) const noexcept {
      return data[((n*n3 + k)*n2 + j)*n1 + i];
    }
    inline  T& operator()(int n, int k, int j, int i)  noexcept {
      return data[((n*n3 + k)*n2 + j)*n1 + i];
    }
  };

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
  
};
#endif
