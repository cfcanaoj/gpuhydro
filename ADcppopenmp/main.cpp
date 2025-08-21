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

using namespace hydro_arrays_mod;

namespace resolution_mod {
  int stepmax{100}; // max step 
  int stepsnap = stepmax/100;
  double time = 0.0e0;
  double dt;
  double timemax = 5.0e0;
  double dtout = 5.0e0/600;
  
  int nx{64}; //! resolution for x
  int ny{64}; //! resolution for y
  int nz{64}; //! resolution for z
  int ngh{2};  //! numeber of ghost mesh
  int itot = nx + 2 * ngh; // 
  int jtot = ny + 2 * ngh; // 
  int ktot = nz + 2 * ngh; // 
  int is = ngh; // 
  int js = ngh; // 
  int ks = ngh; // 
  int ie = is + nx-1;
  int je = js + ny-1;
  int ke = ks + nz-1;
  double xmin(-0.5),xmax(+0.5);
  double ymin(-0.5),ymax(+0.5);
  double zmin(-0.5),zmax(+0.5);
  double dx = (xmax-xmin)/nx;
  double dy = (ymax-ymin)/ny;
  double dz = (zmax-zmin)/nz;
};

namespace hydflux_mod {
#pragma omp declare target
  int mconsv{5}; //!
  HydroArrays<double> U; //! U(mconsv,ktot,jtot,itot)
  int mden{0},mrvx{1},mrvy{2},mrvz{3},meto{4};
  HydroArrays<double> fluxx,fluxy,fluxz;
  int nprim{5}; //!
  HydroArrays<double> P; //! P(nprim,ktot,jtot,itot)
  int nden{0},nvex{1},nvey{2},nvez{3},nene{4};
#pragma omp end declare target 

};

static void SetBoundaryCondition() {
  using namespace resolution_mod;
  using namespace hydflux_mod;

  // x-direction
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je;j++)
	for (int i=1; i<=ngh; i++) {   
	  P(n,k,j,is-i) = P(n,k,j,ie+1-i);
	  P(n,k,j,ie+i) = P(n,k,j,is-1+i);
	  
  }
  

  // y-direction
  for (int n=0; n<nprim; n++)
    for (int k=ks; k<=ke; k++)
      for (int j=1; j<=ngh;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,k,js-j,i) = P(n,k,je+1-j,i);
	  P(n,k,je+j,i) = P(n,k,js-1+j,i);
  }
  
  // z-direction
  for (int n=0; n<nprim; n++)
    for (int k=1; k<=ngh; k++)
      for (int j=js; j<=je;j++)
	for (int i=is; i<=ie; i++) {
	  P(n,ks-k,j,i) = P(n,ke+1-k,j,i);
	  P(n,ke+k,j,i) = P(n,ks-1+k,j,i);
  }
  
}


static void GetNumericalFlux1(){
  using namespace resolution_mod;
  using namespace hydflux_mod;

#pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie+1; ++i) {
	double qL = P(nden,k,j,i-1);
	double qR = P(nden,k,j,i  );
	double ux = 0.5e0*(P(nvex,k,j,i-1)+P(nvex,k,j,i));
	double q_up = (ux >= 0.0) ? qL : qR;
	fluxx(mden,k,j,i) = ux * q_up;
      }
}


static void GetNumericalFlux2(){
  using namespace resolution_mod;
  using namespace hydflux_mod;

#pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je+1; ++j)
      for (int i=is; i<=ie; ++i) {
	double qL = P(nden,k,j-1,i);
	double qR = P(nden,k,j  ,i);
	double uy = 0.5e0*(P(nvey,k,j-1,i)+P(nvey,k,j,i));
	double q_up = (uy >= 0.0) ? qL : qR;
	fluxy(mden,k,j,i) = uy * q_up;
      }
}

static void GetNumericalFlux3(){
  using namespace resolution_mod;
  using namespace hydflux_mod;

#pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke+1; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
	double qL = P(nden,k-1,j,i);
	double qR = P(nden,k  ,j,i);
	double uz = 0.5e0*(P(nvez,k-1,j,i)+P(nvez,k,j,i));
	double q_up = (uz >= 0.0) ? qL : qR;
	fluxz(mden,k,j,i) = uz * q_up;
      }
}

static void UpdateConservU(){
  using namespace resolution_mod;
  using namespace hydflux_mod;

#pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
        double dqx = (fluxx(mden,k,j,i+1) - fluxx(mden,k,j,i)) / dx;
        double dqy = (fluxy(mden,k,j+1,i) - fluxy(mden,k,j,i)) / dy;
        double dqz = (fluxz(mden,k+1,j,i) - fluxz(mden,k,j,i)) / dz;
        U(mden,k,j,i) -= dt * (dqx + dqy + dqz);
      }
}


static void UpdatePrimitvP(){
  using namespace resolution_mod;
  using namespace hydflux_mod;

#pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
         P(nden,k,j,i) = U(mden,k,j,i);
	 //P(nvex,k,j,i) = U(mrvx,k,j,i)/U(mden,k,j,i);
         //P(nvey,k,j,i) = U(mrvy,k,j,i)/U(mden,k,j,i);
         //P(nvez,k,j,i) = U(mrvz,k,j,i)/U(mden,k,j,i);
	 double ekin = 0.5e0*( U(mrvx,k,j,i)*U(mrvx,k,j,i)
			      +U(mrvy,k,j,i)*U(mrvy,k,j,i)
			      +U(mrvz,k,j,i)*U(mrvz,k,j,i))/U(mden,k,j,i);
         P(nene,k,j,i) = U(meto,k,j,i)-ekin;
      }
}

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

static void ControlTimestep(){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  double dtmin;
  const double eps = 1.0e-10;
  dtmin = 1.0e10;
  for (int k=ks; k<=ke; k++)
    for (int j=js; j<=je; j++)
      for (int i=is; i<=ie; i++) {
	double dtminloc = std::min({dx/(P(nvex,k,j,i)+eps)
			           ,dy/(P(nvey,k,j,i)+eps)
				   ,dz/(P(nvez,k,j,i)+eps)});
	dtmin = std::min(dtminloc,dtmin);
      }
  dt = 0.5*dtmin;
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
  using namespace hydflux_mod;
      U.allocate(mconsv,ktot,jtot,itot);
  fluxx.allocate(mconsv,ktot,jtot,itot);
  fluxy.allocate(mconsv,ktot,jtot,itot);
  fluxz.allocate(mconsv,ktot,jtot,itot);
  
      P.allocate(nprim ,ktot,jtot,itot);

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
