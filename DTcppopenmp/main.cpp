/**
 * @file main.cpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <chrono>

#include "hydro_arrays.hpp"
#include "resolution.hpp"
#include "hydro.hpp"
#include "boundary.hpp"

#include "main.hpp"


using namespace hydro_arrays_mod;

static void GenerateProblem(Array4D<double>& P,Array4D<double>& U) {
  using namespace resolution_mod;
  using namespace hydflux_mod;

  double pi;
  double Ahl=0.5e0,Bhl=0.5e0,Chl=0.5e0;
  double k_ini=2.0e0;
  double ekin = 2.0e0;
  double emag = 2.0e0;
  double eint = 1.0e0;
  double d0 = 1.0e0;
  double v0;
  double b0;
  double p0;
  double eps = 1.0e-1;
  double deltax = 0.1e0,deltay = 0.2e0,deltaz = 0.3e0; // randam phase

  pi = acos(-1.0e0);
  
  v0 = sqrt(ekin*2.0e0/d0);
  b0 = sqrt(emag*2.0);
  csiso = sqrt(eint/d0);
  p0 = d0*csiso*csiso;
  
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
	double x = (i-is+0.5e0)*dx;
	double y = (j-js+0.5e0)*dy;
	double z = (k-ks+0.5e0)*dz;
	  
	P(nden,k,j,i) = d0;
	P(nve1,k,j,i) = v0*(  Ahl*sin(2.0*pi*(k_ini*z/(zmax-zmin)+deltaz))
			     +Chl*cos(2.0*pi*(k_ini*y/(ymax-ymin)+deltay)) );
	P(nve2,k,j,i) = v0*(  Bhl*sin(2.0*pi*(k_ini*x/(xmax-xmin)+deltax))
			     +Ahl*cos(2.0*pi*(k_ini*z/(zmax-zmin)+deltaz)) );
	P(nve3,k,j,i) = v0*(  Chl*sin(2.0*pi*(k_ini*y/(ymax-ymin)+deltay))
			     +Bhl*cos(2.0*pi*(k_ini*x/(xmax-xmin)+deltax)) );
	
	P(npre,k,j,i)  =p0;

	P(nbm1,k,j,i) = b0*(  Ahl*sin(2.0*pi*(k_ini*z/(zmax-zmin)+deltaz))
			     +Chl*cos(2.0*pi*(k_ini*y/(ymax-ymin)+deltay)) );
	P(nbm2,k,j,i) = b0*(  Bhl*sin(2.0*pi*(k_ini*x/(xmax-xmin)+deltax))
			     +Ahl*cos(2.0*pi*(k_ini*z/(zmax-zmin)+deltaz)) );
	P(nbm3,k,j,i) = b0*(  Chl*sin(2.0*pi*(k_ini*y/(ymax-ymin)+deltay))
			     +Bhl*cos(2.0*pi*(k_ini*x/(xmax-xmin)+deltax)) );
    };

  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
	P(nden,k,j,i) = P(npre,k,j,i);
	P(ncsp,k,j,i) = csiso;
    };


  
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
    U(mden,k,j,i) = P(nden,k,j,i);
    U(mrv1,k,j,i) = P(nden,k,j,i)*P(nve1,k,j,i);
    U(mrv2,k,j,i) = P(nden,k,j,i)*P(nve2,k,j,i);
    U(mrv3,k,j,i) = P(nden,k,j,i)*P(nve3,k,j,i);
    double ekin = 0.5*P(nden,k,j,i)*( P(nve1,k,j,i)*P(nve1,k,j,i)
	        		     +P(nve2,k,j,i)*P(nve2,k,j,i)
				     +P(nve3,k,j,i)*P(nve3,k,j,i));
    double emag = 0.5              *( P(nbm1,k,j,i)*P(nbm1,k,j,i)
	        		     +P(nbm2,k,j,i)*P(nbm2,k,j,i)
				     +P(nbm3,k,j,i)*P(nbm3,k,j,i));
     U(meto,k,j,i) = P(nene,k,j,i) + ekin + emag;
  };
  
}

void Output1D(){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  static int index = 0;
  FILE *ofile;
  int i;
  char outfile[20];
  int          ret;

  int jc = int((js+je)/2);
  int kc = int((ks+ke)/2);
  ret=sprintf(outfile,"snap/t%05d.dat",index);
  (void)system("mkdir -p snap");
  ofile = fopen(outfile,"w");
  fprintf(ofile,  "# %12.7e\n",time_sim);
  fprintf(ofile,  "# %12s %12s %12s\n","x[cm]","rho[g/cm^3]","v_x[cm/s]");
  for(i=is;i<=ie;i++){
    fprintf(ofile,"  %12.5e %12.5e %12.5e \n",xmin+(i-is+0.5)*dx,P(nden,kc,jc,i),P(nve1,kc,jc,i));
  }
  fclose(ofile);
  index += 1;
}

int main() {
  
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace boundary_mod;
  
  printf("setup grids and fields\n");
  
  AllocateHydroVariables(U,Fx,Fy,Fz,P);
 
  AllocateBoundaryVariables(Xs,Xe,Ys,Ye,Zs,Ze);
  
  printf("grid size for x y z = %i %i %i\n",nx,ny,nz);
  
  GenerateProblem(P,U);
  printf("entering main loop\n");
  int step = 0;
  auto time_begin = std::chrono::high_resolution_clock::now();

  for (step=0;step<stepmax;step++){

    ControlTimestep(); 
    SetBoundaryCondition(P,Xs,Xe,Ys,Ye,Zs,Ze);
    EvaluateCh();
    GetNumericalFlux1(P,Fx);
    GetNumericalFlux2(P,Fy);
    GetNumericalFlux3(P,Fz);
    UpdateConservU(Fx,Fy,Fz,U);
    UpdatePrimitvP(U,P);

    time_sim += dt;
    //printf("dt=%e\n",dt);
    if (step % stepsnap == 0) {
#pragma omp target update from (P.data[0:P.size])
      Output1D();
    }
  }

  //DeallocateHydroVariables(U,Fx,Fy,Fz,P);
  //DeallocateBoundaryVariables(Xs,Xe,Ys,Ye,Zs,Ze);
  
  auto time_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = time_end - time_begin;
  printf("sim time [s]: %e\n", elapsed.count());
  printf("time/count/cell : %e\n", elapsed.count()/(nx*ny*nz)/stepmax);
    
  printf("program has been finished\n");
  return 0;
}
