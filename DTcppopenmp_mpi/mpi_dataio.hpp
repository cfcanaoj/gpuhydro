#ifndef MPI_DATAIO_HPP_
#define MPI_DATAIO_HPP_

#include "mpi_config.hpp"
namespace mpi = mpi_config_mod;
  
using index_t = int;
using size_t  = std::size_t;

namespace mpi_dataio_mod {

  template <typename T>
  class GridArray {
  public:
    T* data = nullptr;
    int nv = 0, n1 = 0;
    size_t size = 0;
    GridArray() = default;
    GridArray(int nv_, int n1_) {
      allocate(nv_,n1_);
    }
    void allocate(int nv_, int n1_) {
      nv = nv_; n1 = n1_;
      size = static_cast<size_t>(nv)*n1;
      data = static_cast<T*>(malloc(sizeof(T) * size));
    }

    inline       T& operator()(int n,int i)        noexcept { return data[(n)*n1+i]; }
    inline const T& operator()(int n,int i) const  noexcept { return data[(n)*n1+i]; }
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

    inline const T& operator()(int n, int k, int j, int i) const noexcept {
      return data[((n*n3 + k)*n2 + j)*n1 + i];
    }
    inline  T& operator()(int n, int k, int j, int i)  noexcept {
      return data[((n*n3 + k)*n2 + j)*n1 + i];
    }
  };

  extern char id[10];
  extern char datadir[10];
  
  extern int ntotal[mpi::ndim];
  extern int npart[mpi::ndim];
  extern int nvars, nvarg;

  extern GridArray<double> gridXout;
  extern GridArray<double> gridYout;
  extern GridArray<double> gridZout;

  extern FieldArray<double> Fieldout;

  void MPIOutputBindary(int timeid);  
}

#endif 
