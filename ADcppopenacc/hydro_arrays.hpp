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

//----------------------------------------------
// 4D 配列: U(n,k,j,i) 変数×3D（SoA 的）
//----------------------------------------------
template <typename T>
class Array4D {
public:
  std::vector<T> data;
  index_t n1{0}, n2{0}, n3{0}, nv{0};
  Array4D() = default;
  Array4D(index_t nv_, index_t n3_, index_t n2_, index_t n1_)
  { resize(nv_, n3_, n2_, n1_); }
  
  void allocate(index_t nv_,index_t n3_, index_t n2_, index_t n1_) {
    nv=nv_;n3 = n3_; n2 = n2_; n1 = n1_;
    data.assign(static_cast<size_t>(nv)*n3*n2*n1, T{});
  }
  
  inline       T& operator()(index_t n, index_t k, index_t j, index_t i)       noexcept {
    #ifndef NDEBUG
      assert(0<=i && i<n1 && 0<=j && j<n2 && 0<=k && k<n3 && 0<=n && n<nv);
    #endif
    // i fastest -> (((n)*n3 + k)*n2 + j)*n1 + i
    return data[ (((static_cast<size_t>(n)*n3 + k)*n2) + j)*n1 + i ];
  }
  inline const T& operator()(index_t n, index_t k, index_t j, index_t i) const noexcept {
    #ifndef NDEBUG
      assert(0<=i && i<n1 && 0<=j && j<n2 && 0<=k && k<n3 && 0<=n && n<nv);
    #endif
    return data[ (((static_cast<size_t>(n)*n3 + k)*n2) + j)*n1 + i ];
  }

  // 次元
  size_t  size() const noexcept { return n1*n2*n3*nv; }

  void fill(const T& v) { std::fill(data.begin(), data.end(), v); }

};

} // namespace hydromod

#endif 
