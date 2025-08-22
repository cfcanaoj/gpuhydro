/**
 * @file main.cpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/
#include <cmath>
#include <cstdio>
#include <algorithm>
#include "hydro_arrays.hpp"
#include "resolution.hpp"
#include "hydro.hpp"
#include "boundary.hpp"

static void GenerateProblem() {
  using namespace resolution_mod;
  using namespace hydflux_mod;

  double cx = 0.5*(xmax-xmin), cy = 0.5*(ymax-ymin), cz = 0.5*(zmax-zmin);
  double sig2 = 0.01; // 分散（適当に）

  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
    double x = (i - is + 0.5) * dx;
    double y = (j - js + 0.5) * dy;
    double z = (k - ks + 0.5) * dz;
    double dx_ = x - cx, dy_ = y - cy, dz_ = z - cz;
    double r2 = dx_*dx_ + dy_*dy_ + dz_*dz_;
    P(nden,k,j,i) = std::exp(-r2 / (2.0*sig2));
  };
  
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
    P(nvex,k,j,i) = 1.0e0;
    P(nvey,k,j,i) = 0.0e0;
    P(nvez,k,j,i) = 0.0e0;
    P(nene,k,j,i) = 0.0e0;    
  };

  
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
    U(mden,k,j,i) = P(nden,k,j,i);
    U(mrvx,k,j,i) = P(nden,k,j,i)*P(nvex,k,j,i);
    U(mrvx,k,j,i) = P(nden,k,j,i)*P(nvey,k,j,i);
    U(mrvx,k,j,i) = P(nden,k,j,i)*P(nvez,k,j,i);
    double ekin = 0.5*P(nden,k,j,i)*( P(nvex,k,j,i)*P(nvex,k,j,i)
	        		     +P(nvey,k,j,i)*P(nvey,k,j,i)
				     +P(nvez,k,j,i)*P(nvez,k,j,i));
     U(meto,k,j,i) = P(nene,k,j,i)+ekin;
  };
#pragma omp target update to (U,P)
  
}

void Output1D(){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  static int index = 0;
  FILE *ofile;
  int i;
  char outfile[20];
  int          ret;

  int j,k;
  j=js;
  k=ks;
  ret=sprintf(outfile,"snap/t%05d.dat",index);
  ofile = fopen(outfile,"w");
  fprintf(ofile,  "# %12.7e\n",time);
  fprintf(ofile,  "# %12s %12s %12s\n","x[cm]","rho[g/cm^3]","v_x[cm/s]");
  for(i=is;i<=ie;i++){
    fprintf(ofile,"  %12.5e %12.5e %12.5e \n",i*dx,P(nden,k,j,i),P(nvex,k,j,i));
  }
  fclose(ofile);
  index += 1;
}

int main() {
  using namespace resolution_mod;

  AllocateVariables();
  // 初期条件
  GenerateProblem();

  int step = 0;

  for (step=0;step<stepmax;step++){

    ControlTimestep(); 
    SetBoundaryCondition();
    GetNumericalFlux1();
    GetNumericalFlux2();
    GetNumericalFlux3();
    UpdateConservU();
    UpdatePrimitvP();

    time += dt;

    if (step % stepsnap == 0) {
      Output1D();
    }
  }

  return 0;
}
