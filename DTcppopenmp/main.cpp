/**
 * @file main.cpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cerrno>

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
  //v0 = 0.0e0;
  //b0 = 0.0e0;
  csiso = sqrt(eint/d0);
  p0 = d0*csiso*csiso;  
  chg = 0.0e0; // this value is not used
#pragma omp target update to ( csiso,chg)

  //printf("cs=%e",csiso);
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
	double x = xmin+(i-is+0.5e0)*dx;
	double y = ymin+(j-js+0.5e0)*dy;
	double z = zmin+(k-ks+0.5e0)*dz;
	  
	P(nden,k,j,i) = d0;
	P(nve1,k,j,i) = v0*(  Ahl*sin(2.0*pi*(k_ini*z/(zmax-zmin)+deltaz))
			     +Chl*cos(2.0*pi*(k_ini*y/(ymax-ymin)+deltay)) );
	P(nve2,k,j,i) = v0*(  Bhl*sin(2.0*pi*(k_ini*x/(xmax-xmin)+deltax))
			     +Ahl*cos(2.0*pi*(k_ini*z/(zmax-zmin)+deltaz)) );
	P(nve3,k,j,i) = v0*(  Chl*sin(2.0*pi*(k_ini*y/(ymax-ymin)+deltay))
			     +Bhl*cos(2.0*pi*(k_ini*x/(xmax-xmin)+deltax)) );
	P(nene,k,j,i)  = eint/d0; //specific internel energy	
	P(npre,k,j,i)  =p0;
	P(ncsp,k,j,i) = csiso;
	
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
     U(meto,k,j,i) = P(nene,k,j,i)*P(nden,k,j,i) + ekin + emag;
  };
#pragma omp target update to ( U.data[0: U.size])
#pragma omp target update to ( P.data[0: P.size])

}

static void GenerateProblem2(Array4D<double>& P,Array4D<double>& U) {
  using namespace resolution_mod;
  using namespace hydflux_mod;

  double eint = 1.0e0;
  double denc = 1.0e0;
  csiso = sqrt(eint/denc);
  chg = 0.0e0;
#pragma omp target update to ( csiso,chg)

  double pres = denc=csiso*csiso;

  double xc = 0.5*(xmax + xmin);
  double yc = 0.5*(ymax + ymin);
  double zc = 0.5*(zmax + zmin);
  double sigma2x = pow((0.1e0*(xmax-xmin)),2);
  double sigma2y = pow((0.1e0*(ymax-ymin)),2);
  double sigma2z = pow((0.1e0*(zmax-zmin)),2);
  
  //printf("cs=%e",csiso);
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
	double x = xmin+(i-is+0.5e0)*dx;
	double y = ymin+(j-js+0.5e0)*dy;
	double z = zmin+(k-ks+0.5e0)*dz;
	  
	P(nden,k,j,i) = denc;
	P(nve1,k,j,i) = 0.3e0;
	P(nve2,k,j,i) = 0.3e0;
	P(nve3,k,j,i) = 0.3e0;
	P(nene,k,j,i)  = 1.0e6*exp(-( pow(x-xc,2)/sigma2x
				     +pow(y-yc,2)/sigma2y
	                             +pow(z-zc,2)/sigma2z )); //specific internel energy	
	P(npre,k,j,i)  =pres;
	P(ncsp,k,j,i) = csiso;
	
	P(nbm1,k,j,i) = 0.0e0;
	P(nbm2,k,j,i) = 0.0e0;
	P(nbm3,k,j,i) = 0.0e0;
	P(nbps,k,j,i) = 0.0e0;
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
     U(meto,k,j,i) = P(nene,k,j,i)*P(nden,k,j,i) + ekin + emag;
  };
#pragma omp target update to ( U.data[0: U.size])
#pragma omp target update to ( P.data[0: P.size])

}

void Output(bool& forcedamp){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  static int index = 0;
  const int gs = 1;
  const int nvar = 9;
  
  static Array4D<double> hydout;
  static bool is_inited = false;


  if(!forcedamp && time_sim < time_out + dtout) return;
  

#pragma omp target update from (P.data[0:P.size])
  
  if (! is_inited){
    (void)system("mkdir -p bindata");
    hydout.allocate(nvar,nz+2*gs,ny+2*gs,nx*2*gs);
    is_inited = true;
  }
  // ---- output text (unf%05d.dat) ----
  char fname_unf[256];
  std::snprintf(fname_unf, sizeof(fname_unf), "bindata/unf%05d.dat", index);
  FILE* fp_unf = std::fopen(fname_unf, "w");
  if (!fp_unf){
    std::fprintf(stderr, "open failed: %s : %s\n", fname_unf, std::strerror(errno));
  }
  
  std::fprintf(fp_unf, "# %.17g %.17g\n", time_sim,dt);
  std::fprintf(fp_unf, "# %d %d\n", nx, gs);
  std::fprintf(fp_unf, "# %d %d\n", ny, gs);
  std::fprintf(fp_unf, "# %d %d\n", nz, gs);
  std::fclose(fp_unf);
  
  // ---- output data (bin%05d.dat) ----

  for (int k=ks-gs;k<=ke+gs;k++)
    for (int j=js-gs;j<=je+gs;j++)
      for (int i=is-gs;i<=ie+gs;i++){
	hydout(0,k-1,j-1,i-1) = P(nden,k,j,i);
	hydout(1,k-1,j-1,i-1) = P(nve1,k,j,i);
	hydout(2,k-1,j-1,i-1) = P(nve2,k,j,i);
	hydout(3,k-1,j-1,i-1) = P(nve3,k,j,i);
	hydout(4,k-1,j-1,i-1) = P(nbm1,k,j,i);
	hydout(5,k-1,j-1,i-1) = P(nbm2,k,j,i);
	hydout(6,k-1,j-1,i-1) = P(nbm3,k,j,i);
	hydout(7,k-1,j-1,i-1) = P(nbps,k,j,i);
	hydout(8,k-1,j-1,i-1) = P(npre,k,j,i); 
  }
  
  char fname_bin[256];
  std::snprintf(fname_bin, sizeof(fname_bin), "bindata/bin%05d.dat",index);
  FILE* fp_bin = std::fopen(fname_bin, "wb");

  auto wbin = [&](const void* ptr, size_t n) {
		const size_t wrote = std::fwrite(ptr, sizeof(double), n, fp_bin);
		if (wrote != n) {
		  std::fclose(fp_bin);
		  if (!fp_bin){
		    std::fprintf(stderr, "open failed: %s : %s\n", fname_bin, std::strerror(errno));
		  }
		}
	      };

  wbin(hydout.data, hydout.size);
  std::fclose(fp_bin);

  index += 1;
}


void Output1D(bool& forcedamp){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  static int index = 0;  
  static bool is_inited = false;

  if(!forcedamp && time_sim < time_out + dtout) return;

  printf("output index=%i, time=%e",index,time_sim);

#pragma omp target update from (P.data[0:P.size])
  int ic,jc,kc;

  ic = int((is+ie)/2);
  jc = int((js+je)/2);
  kc = int((ks+ke)/2);
  
  if (! is_inited){
    (void)system("mkdir -p snap");
    is_inited = true;
  }
  // for x
  char fname[256];
  std::snprintf(fname, sizeof(fname), "snap/snx%05d.dat", index);
  FILE* fp = std::fopen(fname, "w");
  if (!fp){
    std::fprintf(stderr, "open failed: %s : %s\n", fname, std::strerror(errno));
  }
  for (int i=is-2;i<=ie+2;i++){
    std::fprintf(fp, "%e %e %e \n", xmin+(i-is+0.5)*dx,P(nden,kc,jc,i),P(nene,kc,jc,i) );
  }
  std::fclose(fp);
  
  // for y
  std::snprintf(fname, sizeof(fname), "snap/sny%05d.dat", index);
  fp = std::fopen(fname, "w");
  if (!fp){
    std::fprintf(stderr, "open failed: %s : %s\n", fname, std::strerror(errno));
  }
  for (int j=js-2;j<=je+2;j++){
    std::fprintf(fp, "%e %e %e \n", ymin+(j-js+0.5)*dy,P(nden,kc,j,ic),P(nene,kc,j,ic) );
  }
  std::fclose(fp);

  // for z
  std::snprintf(fname, sizeof(fname), "snap/snz%05d.dat", index);
  fp = std::fopen(fname, "w");
  if (!fp){
    std::fprintf(stderr, "open failed: %s : %s\n", fname, std::strerror(errno));
  }
  for (int k=ks-2;k<=ke+2;k++){
    std::fprintf(fp, "%e %e %e \n", zmin+(k-ks+0.5)*dz,P(nden,k,jc,ic),P(nene,k,jc,ic) );
  }
  std::fclose(fp);



  
  index += 1;
}


int main() {
  
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace boundary_mod;
  const bool NoOutput = true;
  static bool is_final = false;
  printf("setup grids and fields\n");
  
  AllocateHydroVariables(U,Fx,Fy,Fz,P);
 
  AllocateBoundaryVariables(Xs,Xe,Ys,Ye,Zs,Ze);
  
  printf("grid size for x y z = %i %i %i\n",nx,ny,nz);
  
  GenerateProblem(P,U);
  //GenerateProblem2(P,U);
  printf("entering main loop\n");
  int step = 0;
  auto time_beg = std::chrono::high_resolution_clock::now();

  for (step=0;step<stepmax;step++){

    ControlTimestep(); 
    if (step%300 ==0 && !NoOutput) printf("step=%i time=%e dt=%e\n",step,time_sim,dt);
    //printf("step=%i time=%e dt=%e\n",step,time_sim,dt);
    SetBoundaryCondition(P,Xs,Xe,Ys,Ye,Zs,Ze);
    EvaluateCh();
    GetNumericalFlux1(P,Fx);
    GetNumericalFlux2(P,Fy);
    GetNumericalFlux3(P,Fz);
    UpdateConservU(Fx,Fy,Fz,U);
    DampPsi(U);
    UpdatePrimitvP(U,P);

    time_sim += dt;
    //printf("dt=%e\n",dt);
    if (! NoOutput) Output(is_final);
    //if (! NoOutput) Output1D(is_final);

    if(time_sim > time_max) break;
    
  }

  //DeallocateHydroVariables(U,Fx,Fy,Fz,P);
  //DeallocateBoundaryVariables(Xs,Xe,Ys,Ye,Zs,Ze);
  
  auto time_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = time_end - time_beg;
  printf("exiting main loop time=%e, step=%i\n",time_sim,step);
  printf("sim time [s]: %e\n", elapsed.count());
  printf("time/count/cell : %e\n", elapsed.count()/(nx*ny*nz)/stepmax);

  is_final = true;
  Output(is_final);
  //Output1D(is_final);
  
  printf("program has been finished\n");
  return 0;
}
