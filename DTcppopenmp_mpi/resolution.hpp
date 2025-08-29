/**
 * @file resolution.hpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/

namespace resolution_mod {
  inline constexpr int stepmax{3000}; // max step 
  inline constexpr int stepsnap = stepmax/100;
  inline double time_sim = 0.0e0;
  inline double time_out = 0.0e0;
  inline double dt;
  inline constexpr double time_max = 3.0e0;
  inline constexpr double dtout = time_max/100;
  
  inline constexpr int ngrid1{128}; //! resolution for x
  inline constexpr int ngrid2{64}; //! resolution for y
  inline constexpr int ngrid3{128}; //! resolution for z
  inline constexpr int ngh{2};  //! number of ghost mesh
  inline constexpr int itot = ngrid1 + 2 * ngh+1; //! Like ZEUS-2D
  inline constexpr int jtot = ngrid2 + 2 * ngh+1; // 
  inline constexpr int ktot = ngrid3 + 2 * ngh+1; // 
  inline constexpr int is = ngh; //! |0 1 | 2 for ngh =2
  inline constexpr int js = ngh; // 
  inline constexpr int ks = ngh; // 
  inline constexpr int ie = is + ngrid1-1; //! 65 | 66 67 | for nx =64 
  inline constexpr int je = js + ngrid2-1;
  inline constexpr int ke = ks + ngrid3-1;
  inline constexpr double x1min(-0.5),x1max(+0.5);
  inline constexpr double x2min(-0.5),x2max(+0.5);
  inline constexpr double x3min(-0.5),x3max(+0.5);
};
