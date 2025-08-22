/**
 * @file resolution.hpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/

namespace resolution_mod {
  inline constexpr int stepmax{100}; // max step 
  inline constexpr int stepsnap = stepmax/100;
  inline double time = 0.0e0;
  inline double dt;
  inline constexpr double timemax = 5.0e0;
  inline constexpr double dtout = 5.0e0/600;
  
  inline constexpr int nx{64}; //! resolution for x
  inline constexpr int ny{64}; //! resolution for y
  inline constexpr int nz{64}; //! resolution for z
  inline constexpr int ngh{2};  //! numeber of ghost mesh
  inline constexpr int itot = nx + 2 * ngh; // 
  inline constexpr int jtot = ny + 2 * ngh; // 
  inline constexpr int ktot = nz + 2 * ngh; // 
  inline constexpr int is = ngh; // 
  inline constexpr int js = ngh; // 
  inline constexpr int ks = ngh; // 
  inline constexpr int ie = is + nx-1;
  inline constexpr int je = js + ny-1;
  inline constexpr int ke = ks + nz-1;
  inline constexpr double xmin(-0.5),xmax(+0.5);
  inline constexpr double ymin(-0.5),ymax(+0.5);
  inline constexpr double zmin(-0.5),zmax(+0.5);
  inline constexpr double dx = (xmax-xmin)/nx;
  inline constexpr double dy = (ymax-ymin)/ny;
  inline constexpr double dz = (zmax-zmin)/nz;
};
