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
// 3D 配列: A(k,j,i) で i が連続（SoA のスカラー場）
//----------------------------------------------
template <typename T>
class Array3D {
public:
  Array3D() = default;

  Array3D(index_t n3, index_t n2, index_t n1)
  { resize(n3, n2, n1); }

  void resize(index_t n3, index_t n2, index_t n1) {
    n3_ = n3; n2_ = n2; n1_ = n1;
    data_.assign(static_cast<size_t>(n3_)*n2_*n1_, T{});
  }

  // 連続メモリ先頭ポインタ
  T* data() noexcept { return data_.data(); }
  const T* data() const noexcept { return data_.data(); }

  // 要素アクセス (境界チェックはデバッグ時のみ)
  inline       T& operator()(index_t k, index_t j, index_t i)       noexcept {
    #ifndef NDEBUG
      assert(0<=i && i<n1_ && 0<=j && j<n2_ && 0<=k && k<n3_);
    #endif
    return data_[ (static_cast<size_t>(k)*n2_ + j)*n1_ + i ];
  }
  inline const T& operator()(index_t k, index_t j, index_t i) const noexcept {
    #ifndef NDEBUG
      assert(0<=i && i<n1_ && 0<=j && j<n2_ && 0<=k && k<n3_);
    #endif
    return data_[ (static_cast<size_t>(k)*n2_ + j)*n1_ + i ];
  }

  // 次元
  index_t n1() const noexcept { return n1_; } // i
  index_t n2() const noexcept { return n2_; } // j
  index_t n3() const noexcept { return n3_; } // k
  size_t size() const noexcept { return data_.size(); }

  void fill(const T& v) { std::fill(data_.begin(), data_.end(), v); }

private:
  index_t n1_{0}, n2_{0}, n3_{0}; // i fastest, then j, then k
  std::vector<T> data_;
};

//----------------------------------------------
// 4D 配列: U(n,k,j,i) 変数×3D（SoA 的）
//----------------------------------------------
template <typename T>
class Array4D {
public:
  Array4D() = default;
  Array4D(index_t nv, index_t n3, index_t n2, index_t n1)
  { resize(nv, n3, n2, n1); }

  void resize(index_t nv, index_t n3, index_t n2, index_t n1) {
    nv_ = nv; n3_ = n3; n2_ = n2; n1_ = n1;
    data_.assign(static_cast<size_t>(nv_)*n3_*n2_*n1_, T{});
  }

  inline       T& operator()(index_t n, index_t k, index_t j, index_t i)       noexcept {
    #ifndef NDEBUG
      assert(0<=i && i<n1_ && 0<=j && j<n2_ && 0<=k && k<n3_ && 0<=n && n<nv_);
    #endif
    // i fastest -> (((n)*n3 + k)*n2 + j)*n1 + i
    return data_[ (((static_cast<size_t>(n)*n3_ + k)*n2_) + j)*n1_ + i ];
  }
  inline const T& operator()(index_t n, index_t k, index_t j, index_t i) const noexcept {
    #ifndef NDEBUG
      assert(0<=i && i<n1_ && 0<=j && j<n2_ && 0<=k && k<n3_ && 0<=n && n<nv_);
    #endif
    return data_[ (((static_cast<size_t>(n)*n3_ + k)*n2_) + j)*n1_ + i ];
  }

  // 次元
  index_t nv() const noexcept { return nv_; } // 変数数
  index_t n1() const noexcept { return n1_; } // i
  index_t n2() const noexcept { return n2_; } // j
  index_t n3() const noexcept { return n3_; } // k
  size_t  size() const noexcept { return data_.size(); }

  void fill(const T& v) { std::fill(data_.begin(), data_.end(), v); }

  T* data() noexcept { return data_.data(); }
  const T* data() const noexcept { return data_.data(); }

private:
  index_t nv_{0}, n1_{0}, n2_{0}, n3_{0};
  std::vector<T> data_;
};

//----------------------------------------------
// 便利ラッパ：保存量/原始量（(n,k,j,i)）
//----------------------------------------------
template <typename T>
class HydroArrays {
public:
  void allocate(index_t nvar, index_t nx1, index_t nx2, index_t nx3) {
    const index_t n1 = nx1;
    const index_t n2 = nx2;
    const index_t n3 = nx3;
    U_.resize(nvar, n3, n2, n1);
  }
  inline       T& operator()(index_t n, index_t k, index_t j, index_t i)       noexcept {
      return U_(n,k,j,i);
  }
  inline const T& operator()(index_t n, index_t k, index_t j, index_t i) const noexcept {
      return U_(n,k,j,i);
  }
  Array4D<T>&       data()       noexcept { return U_; }
  const Array4D<T>& data() const noexcept { return U_; }
  index_t n1() const noexcept {return U_.n1();}
  index_t n2() const noexcept {return U_.n2();}
  index_t n3() const noexcept {return U_.n3();}
  index_t nv() const noexcept {return U_.nv();}
  index_t size() const noexcept {return U_.n1()*U_.n2()*U_.n3()*U_.nv();}
private:
  Array4D<T> U_;
};


} // namespace hydromod

#endif 
