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
// 基本：アクティブ領域/格子メタ
//----------------------------------------------
struct MeshSize {
  index_t nx1{0}, nx2{0}, nx3{0}; // アクティブ(計算領域)の各次元
  index_t nghost{1};              // ゴースト幅（全方向同一と仮定）
};

struct ActiveZone {
  index_t is{0}, ie{-1};
  index_t js{0}, je{-1};
  index_t ks{0}, ke{-1};
};

// ヘルパ：MeshSize -> 確保サイズ & ActiveZone
inline std::pair<ActiveZone, MeshSize> make_layout_with_ghosts(index_t nx1, index_t nx2, index_t nx3, index_t nghost) {
  MeshSize ms{nx1, nx2, nx3, nghost};
  ActiveZone az;
  az.is =          nghost    ;     az.ie =          nghost + nx1 - 1     ;
  az.js = (nx2>0)? nghost : 0;     az.je = (nx2>0)? nghost + nx2 - 1 : -1;
  az.ks = (nx3>0)? nghost : 0;     az.ke = (nx3>0)? nghost + nx3 - 1 : -1;
  return {az, ms};
}

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
    // ゴースト込み確保サイズ
    const index_t n1 = nx1;
    const index_t n2 = nx2;
    const index_t n3 = nx3;
    U_.resize(nvar, n3, n2, n1);
  }

  Array4D<T>&       data()       noexcept { return U_; }
  const Array4D<T>& data() const noexcept { return U_; }

private:
  Array4D<T> U_;
};

//----------------------------------------------
// 軽量ユーティリティ
//----------------------------------------------
template <typename Fn>
inline void for_each_cell(const ActiveZone& az, Fn&& fn) {
  for (index_t k = az.ks; k <= az.ke; ++k)
    for (index_t j = az.js; j <= az.je; ++j)
      for (index_t i = az.is; i <= az.ie; ++i)
        fn(k, j, i);
}

} // namespace hydromod

#endif 
