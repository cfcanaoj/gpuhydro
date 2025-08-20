#include <cmath>
#include <cstdio>
#include <algorithm>
#include "hydro_arrays.hpp"

using namespace hydro_arrays_mod;

// ======================= ユーザー設定 =========================
namespace parameters_mod {
  int nx{128}, ny{128}, nz{128};
  int ngh{2};                 // Ghost Mesh
  double Lx{1.0}, Ly{1.0}, Lz{1.0};
  double ux{0.5}, uy{0.2}, uz{0.1};  // 一様移流速度
  double cfl{0.4};
  double t_end{1.0};
  int output_every{50};       // 何ステップごとに統計出力
  int is = ngh + nx-1;
  int js = ngh + ny-1;
  int ks = ngh + nz-1;
  int ie = is + ngh;
  int je = js + ngh;
  int ke = ks + ngh;
};
// =============================================================

// 周期境界（全方向）: U(n=0,k,j,i) をコピー
static void apply_periodic(HydroArrays<double>& U) {
  using namespace parameters_mod;
  auto& A = U.data();
  const int n1 = nx;
  const int n2 = ny;
  const int n3 = nz;

  // x方向
  for (int k=0; k<n3; ++k)
    for (int j=0; j<n2; ++j) {
      for (int g=0; g<ngh; ++g) {
        A(0,k,j,            g           ) = A(0,k,j, n1 - 2*ngh + g); // 左ハロー
        A(0,k,j, n1 - ngh + g)            = A(0,k,j,      ngh + g ); // 右ハロー
      }
    }
  // y方向（ny==1 のスラブにも対応）
  for (int k=0; k<n3; ++k) {
    for (int i=0; i<n1; ++i) {
      for (int g=0; g<ngh; ++g) {
        int jL = g;
        int jR = n2 - ngh + g;
        int j1 = n2 - 2*ngh + g;
        int j2 = ngh + g;
        A(0,k, jL, i) = A(0,k, j1, i);
        A(0,k, jR, i) = A(0,k, j2, i);
      }
    }
  }
  // z方向（nz==1 のスラブにも対応）
  for (int j=0; j<n2; ++j) {
    for (int i=0; i<n1; ++i) {
      for (int g=0; g<ngh; ++g) {
        int kL = g;
        int kR = n3 - ngh + g;
        int k1 = n3 - 2*ngh + g;
        int k2 = ngh + g;
        A(0, kL, j, i) = A(0, k1, j, i);
        A(0, kR, j, i) = A(0, k2, j, i);
      }
    }
  }
}

// 面心フラックスを一次風上で計算（F = u * q_upwind）
static void compute_face_fluxes(
    const HydroArrays<double>& U,
    HydroArrays<double>& fluxx,
    HydroArrays<double>& fluxy,
    HydroArrays<double>& fluxz,
    double ux, double uy, double uz)
{
  using namespace parameters_mod;
  auto& A = U.data();
  // x1 面：i+1/2 におけるフラックス
  {
    auto& Fx = fluxx.data(); // (k,j,i+1) のサイズ
    printf("p4\n");
#pragma omp target teams distribute parallel for collapse(3)
    for (int k=ks; k<=ke; ++k)
      for (int j=js; j<=je; ++j)
        for (int i=is; i<=ie+1; ++i) {
          // 風上値
          double qL = A(0,k,j,i-1); // 左セル
          double qR = A(0,k,j,i  ); // 右セル
          double q_up = (ux >= 0.0) ? qL : qR;
          Fx(0,k,j,i) = ux * q_up;
        }
  }
    printf("p5\n");
  // x2 面：j+1/2
  {
    auto& Fy = fluxy.data();
#pragma omp target teams distribute parallel for collapse(3)
    for (int k=ks; k<=ke; ++k)
      for (int j=js; j<=je+1; ++j)
        for (int i=is; i<=ie; ++i) {
          double qB = A(0,k,j-1,i); // 下（minus）
          double qT = A(0,k,j  ,i); // 上（plus）
          double q_up = (uy >= 0.0) ? qB : qT;
          Fy(0,k,j,i) = uy * q_up;
        }
  }
  // x3 面：k+1/2
  {
    auto& Fz = fluxz.data();
#pragma omp target teams distribute parallel for collapse(3)
    for (int k=ks; k<=ke+1; ++k)
      for (int j=js; j<=je; ++j)
        for (int i=is; i<=ie; ++i) {
          double qB = A(0,k-1,j,i);
          double qT = A(0,k  ,j,i);
          double q_up = (uz >= 0.0) ? qB : qT;
          Fz(0,k,j,i) = uz * q_up;
        }
  }
}

// セル中心更新（FV：q^{n+1} = q^n - dt/dx * (F_{i+1/2}-F_{i-1/2}) ...）
static void update_cell_center(
    HydroArrays<double>& U,
    const HydroArrays<double>& fluxx,
    const HydroArrays<double>& fluxy,
    const HydroArrays<double>& fluxz,
    double dt, double dx, double dy, double dz)
{
  using namespace parameters_mod;
  auto& A = U.data();

  const auto& Fx = fluxx.data();
  const auto& Fy = fluxy.data();
  const auto& Fz = fluxz.data();

#pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
        double dqx = (Fx(0,k,j,i+1) - Fx(0,k,j,i)) / dx;
        double dqy = (Fy(0,k,j+1,i) - Fy(0,k,j,i)) / dy;
        double dqz = (Fz(0,k+1,j,i) - Fz(0,k,j,i)) / dz;
        A(0,k,j,i) -= dt * (dqx + dqy + dqz);
      }
}

// 初期条件：3D ガウス（周期領域にマッチ）
static void init_gaussian(HydroArrays<double>& U) {
  using namespace parameters_mod;
  auto& A = U.data();

  double cx = 0.5*Lx, cy = 0.5*Ly, cz = 0.5*Lz;
  double sig2 = 0.01; // 分散（適当に）

  double dx = Lx / nx, dy = Ly / ny, dz = Lz / nz;

  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
    double x = (i - is + 0.5) * dx;
    double y = (j - js + 0.5) * dy;
    double z = (k - ks + 0.5) * dz;
    double dx_ = x - cx, dy_ = y - cy, dz_ = z - cz;
    double r2 = dx_*dx_ + dy_*dy_ + dz_*dz_;
    A(0,k,j,i) = std::exp(-r2 / (2.0*sig2));
      };
}

// 統計出力（L1/L2 最大など）
static void print_stats(const HydroArrays<double>& U, const char* tag) {
  using namespace parameters_mod;
  const auto& A = U.data();

  double sum = 0.0, sum2 = 0.0, vmax = 0.0;
  std::size_t ncell = 0;

for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
    double v = A(0,k,j,i);
    sum  += v;
    sum2 += v*v;
    vmax  = std::max(vmax, std::abs(v));
    ++ncell;
  };

  double mean = sum / (double)ncell;
  double rms  = std::sqrt(sum2 / (double)ncell);
  std::printf("[%s] mean=% .6e  rms=% .6e  max|q|=% .6e\n", tag, mean, rms, vmax);
}

int main() {
  using namespace parameters_mod;

  // 格子幅
  double dx = Lx / nx;
  double dy = Ly / ny;
  double dz = Lz / nz;

  // 保存量（nvar=1）と面中心フラックスの確保
  HydroArrays<double> U;
  U.allocate(/*nvar*/1, nx+2*ngh+1,ny+2*ngh+1, nz+2*ngh+1);

  HydroArrays<double> fluxx,fluxy,fluxz;
  fluxx.allocate(/*nvar*/1, nx+2*ngh+1,ny+2*ngh+1, nz+2*ngh+1);
  fluxy.allocate(/*nvar*/1, nx+2*ngh+1,ny+2*ngh+1, nz+2*ngh+1);
  fluxz.allocate(/*nvar*/1, nx+2*ngh+1,ny+2*ngh+1, nz+2*ngh+1);

  // 初期条件
  init_gaussian(U);
  apply_periodic(U);
  print_stats(U, "init");

#pragma omp target enter data map(to:U,fluxx,fluxy,fluxz)
  // 目標時間までステップ
  double t = 0.0;
  int step = 0;
  const double umax =
      std::max({std::abs(ux), std::abs(uy), std::abs(uz), 1e-14});

  while (t < t_end) {
    // CFL に基づく Δt
    double dt_x = (std::abs(ux) > 0) ? cfl * dx / std::abs(ux) : 1e9;
    double dt_y = (std::abs(uy) > 0) ? cfl * dy / std::abs(uy) : 1e9;
    double dt_z = (std::abs(uz) > 0) ? cfl * dz / std::abs(uz) : 1e9;
    double dt = std::min({dt_x, dt_y, dt_z, t_end - t});
    if (dt <= 0) break;

    // 周期境界を更新 → 面心フラックス → セル更新
    printf("p1\n");
    apply_periodic(U);
    printf("p2\n");
    compute_face_fluxes(U, fluxx,fluxy,fluxz, ux, uy, uz);
    printf("p3\n");
    update_cell_center( U, fluxx,fluxy,fluxz, dt, dx, dy, dz);

    t += dt;
    ++step;

    if (step % output_every == 0 || t >= t_end) {
      char buf[64];
      std::snprintf(buf, sizeof(buf), "t=%.4f step=%d", t, step);
      print_stats(U, buf);
    }
  }

#pragma omp target exit data map(from:U,fluxx,fluxy,fluxz)
  // 終了時の簡易チェック（値域）
  print_stats(U, "final");
  return 0;
}
