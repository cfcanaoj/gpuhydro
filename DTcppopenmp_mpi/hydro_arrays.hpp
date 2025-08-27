#ifndef HYDRO_ARRAYS_HPP_
#define HYDRO_ARRAYS_HPP_

// C++ headers

#include <vector>
#include <cassert>
#include <cstddef>
#include <algorithm>
#include <utility>
#include <type_traits>

//
// hydro_arrays.hpp
// - 3D/4D 配列 (i が連続；C 配列的) : (k,j,i) と (n,k,j,i)
// - ゴーストセル込み確保とアクティブ領域インデックス
// - 面中心 3D 配列 (x1/x2/x3 方向で +1 されたサイズを持つ)
// 依存：C++17 以上（std::size_t, std::vector）
//

namespace hydro_arrays_mod {

using index_t = int;     // 必要なら 32/64 を切り替え
using size_t  = std::size_t;
template <typename T>
class Grid3D {
public:
  T* x1a_data = nullptr;
  T* x1b_data = nullptr;
  T* x2a_data = nullptr;
  T* x2b_data = nullptr;
  T* x3a_data = nullptr;
  T* x3b_data = nullptr;
  int n3 = 0, n2 = 0, n1 = 0;
  size_t size = 0;
  Grid3D() = default;
  Grid3D(int n3_, int n2_, int n1_) {
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

  
//----------------------------------------------
// 4D 配列: U(n,k,j,i) 変数×3D（SoA 的）
//----------------------------------------------
template <typename T>
class Array4D {
public:
  T* data = nullptr;
  int nv = 0, n3 = 0, n2 = 0, n1 = 0;
  size_t size = 0;
  Array4D() = default;
  Array4D(int _nv, int _n3, int _n2, int _n1) {
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

} // namespace hydromod

#endif 
