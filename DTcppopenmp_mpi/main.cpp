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

#include "mpi_config.hpp"

#include "mpi_dataio.hpp"

using namespace hydro_arrays_mod;

static void GenerateGrid(GridArray<double>& G) {
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace mpi_config_mod;
  
  double x1minloc,x1maxloc;
  double x2minloc,x2maxloc;
  double x3minloc,x3maxloc;
  double dx1,dx2,dx3;
   
  x1minloc = x1min + (x1max-x1min)/ntiles[dir1]* coords[dir1];
  x1maxloc = x1min + (x1max-x1min)/ntiles[dir1]*(coords[dir1]+1);
  
  dx1 = (x1maxloc-x1minloc)/double(ngrid1);
  for(int i=is-ngh;i<= ie+ngh+1;i++){
    G.x1a(i) = dx1*(i-(ngh+1))+x1minloc;
  }
  for(int i=is-ngh;i<= ie+ngh;i++){
    G.x1b(i) = 0.5e0*(G.x1a(i+1)+G.x1a(i));
  }

  x2minloc = x2min + (x2max-x2min)/ntiles[dir2]* coords[dir2];
  x2maxloc = x2min + (x2max-x2min)/ntiles[dir2]*(coords[dir2]+1);
  dx2=(x2maxloc-x2minloc)/double(ngrid2);
  for(int j=js-ngh;j<= je+ngh+1;j++){
    G.x2a(j) = dx2*(j-(ngh+1))+x2minloc;
  }
  for(int j=js-ngh;j<= je+ngh;j++){
    G.x2b(j) = 0.5e0*(G.x2a(j+1)+G.x2a(j));
  }

  x3minloc = x3min + (x3max-x3min)/ntiles[dir3]* coords[dir3];
  x3maxloc = x3min + (x3max-x3min)/ntiles[dir3]*(coords[dir3]+1);
  dx3=(x3maxloc-x3minloc)/double(ngrid3);
  for(int k=ks-ngh;k<= ke+ngh+1;k++){
    G.x3a(k) = dx3*(k-(ngh+1))+x3minloc;
  }
  for(int k=ks-ngh;k<= ke+ngh;k++){
    G.x3b(k) = 0.5e0*(G.x3a(k+1)+G.x3a(k));
  }

#pragma omp target update to ( G.x1a_data[0:G.n1],G.x1b_data[0:G.n1])
#pragma omp target update to ( G.x2a_data[0:G.n2],G.x2b_data[0:G.n2])
#pragma omp target update to ( G.x3a_data[0:G.n3],G.x3b_data[0:G.n3])

  
}
static void GenerateProblem(GridArray<double>& G,FieldArray<double>& P,FieldArray<double>& U) {
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
	double x = G.x1b(i);
	double y = G.x2b(j);
	double z = G.x3b(k);
	  
	P(nden,k,j,i) = d0;
	P(nve1,k,j,i) = v0*(  Ahl*sin(2.0*pi*(k_ini*z/(x3max-x3min)+deltaz))
			     +Chl*cos(2.0*pi*(k_ini*y/(x2max-x2min)+deltay)) );
	P(nve2,k,j,i) = v0*(  Bhl*sin(2.0*pi*(k_ini*x/(x1max-x1min)+deltax))
			     +Ahl*cos(2.0*pi*(k_ini*z/(x3max-x3min)+deltaz)) );
	P(nve3,k,j,i) = v0*(  Chl*sin(2.0*pi*(k_ini*y/(x2max-x2min)+deltay))
			     +Bhl*cos(2.0*pi*(k_ini*x/(x1max-x1min)+deltax)) );
	P(nene,k,j,i)  = eint/d0; //specific internel energy	
	P(npre,k,j,i)  =p0;
	P(ncsp,k,j,i) = csiso;
	
	P(nbm1,k,j,i) = b0*(  Ahl*sin(2.0*pi*(k_ini*z/(x3max-x3min)+deltaz))
			     +Chl*cos(2.0*pi*(k_ini*y/(x2max-x2min)+deltay)) );
	P(nbm2,k,j,i) = b0*(  Bhl*sin(2.0*pi*(k_ini*x/(x1max-x1min)+deltax))
			     +Ahl*cos(2.0*pi*(k_ini*z/(x3max-x3min)+deltaz)) );
	P(nbm3,k,j,i) = b0*(  Chl*sin(2.0*pi*(k_ini*y/(x2max-x2min)+deltay))
			     +Bhl*cos(2.0*pi*(k_ini*x/(x1max-x1min)+deltax)) );
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

static void GenerateProblem2(GridArray<double>& G,FieldArray<double>& P,FieldArray<double>& U) {
  using namespace resolution_mod;
  using namespace hydflux_mod;

  double eint = 1.0e0;
  double denc = 1.0e0;
  csiso = sqrt(eint/denc);
  chg = 0.0e0;
#pragma omp target update to ( csiso,chg)

  double pres = denc*csiso*csiso;

  double xc = 0.5*(x1max + x1min);
  double yc = 0.5*(x2max + x2min);
  double zc = 0.5*(x3max + x3min);
  double sigma2x = pow((0.1e0*(x1max-x1min)),2);
  double sigma2y = pow((0.1e0*(x2max-x2min)),2);
  double sigma2z = pow((0.1e0*(x3max-x3min)),2);
  
  //printf("cs=%e",csiso);
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
	double x = G.x1b(i);
	double y = G.x2b(j);
	double z = G.x3b(k);
	  
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
     
     U(mbm1,k,j,i) = P(nbm1,k,j,i);
     U(mbm2,k,j,i) = P(nbm2,k,j,i);
     U(mbm3,k,j,i) = P(nbm3,k,j,i);
     U(mbps,k,j,i) = P(nbps,k,j,i);
  };
#pragma omp target update to ( U.data[0: U.size])
#pragma omp target update to ( P.data[0: P.size])

}

void Output1D(bool& forcedamp){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace mpi_config_mod;
  static int index = 0;  
  static bool is_inited = false;
  const int dir1=0, dir2=1, dir3=2;
  if(!forcedamp && time_sim < time_out + dtout) return;

  if(myid_w==0) printf("output index=%i, time=%e \n",index,time_sim);

#pragma omp target update from (P.data[0:P.size])
  int ic,jc,kc;

  // assuming division is up to 2
  if(ntiles[dir1] == 1) {
    ic = int((is+ie)/2);
  }else{
    if(coords[dir1] ==0) ic = ie;
    if(coords[dir1] ==1) ic = is;
  }

  if(ntiles[dir2] == 1) {
    jc = int((js+je)/2);
  }else{
    if(coords[dir2] ==0) jc = je;
    if(coords[dir2] ==1) jc = js;
  }

  
  if(ntiles[dir3] == 1) {
    kc = int((ks+ke)/2);
  }else{
    if(coords[dir3] ==0) kc = ke;
    if(coords[dir3] ==1) kc = ks;
  }
  
  if (! is_inited){
    (void)system("mkdir -p snap");
    is_inited = true;
  }
  // for x
  char fname[256];
  std::snprintf(fname, sizeof(fname), "snap/snx_%02d_%05d.dat",myid_w, index);
  FILE* fp = std::fopen(fname, "w");
  if (!fp){
    std::fprintf(stderr, "open failed: %s : %s\n", fname, std::strerror(errno));
  }
  for (int i=is-2;i<=ie+2;i++){
    std::fprintf(fp, "%e %e %e \n", G.x1b(i),P(nden,kc,jc,i),P(nene,kc,jc,i) );
  }
  std::fclose(fp);
  
  // for y
  std::snprintf(fname, sizeof(fname), "snap/sny_%02d_%05d.dat",myid_w, index);
  fp = std::fopen(fname, "w");
  if (!fp){
    std::fprintf(stderr, "open failed: %s : %s\n", fname, std::strerror(errno));
  }
  for (int j=js-2;j<=je+2;j++){
    std::fprintf(fp, "%e %e %e \n", G.x2b(j),P(nden,kc,j,ic),P(nene,kc,j,ic) );
  }
  std::fclose(fp);

  // for z
  std::snprintf(fname, sizeof(fname), "snap/snz_%02d_%05d.dat",myid_w, index);
  fp = std::fopen(fname, "w");
  if (!fp){
    std::fprintf(stderr, "open failed: %s : %s\n", fname, std::strerror(errno));
  }
  for (int k=ks-2;k<=ke+2;k++){
    std::fprintf(fp, "%e %e %e \n", G.x3b(k),P(nden,k,jc,ic),P(nene,k,jc,ic) );
  }
  std::fclose(fp);
  
  index += 1;
}

void Output(bool& forcedamp){
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace mpi_config_mod;
  namespace io = mpi_dataio_mod;
  static int index = 0;
  
  static bool is_inited = false;

  
  if(!forcedamp && time_sim < time_out + dtout) return;
  

#pragma omp target update from (P.data[0:P.size])
  
  if (! is_inited){
    std::strcpy(io::id,"DT");
    std::strcpy(io::datadir,"./bindata/");
    (void)system("mkdir -p bindata");

    io::ntotal[dir1] = ntiles[dir1]*ngrid1;
    io::ntotal[dir2] = ntiles[dir2]*ngrid2;
    io::ntotal[dir3] = ntiles[dir3]*ngrid3;
    
    io::npart[dir1]  = ngrid1;
    io::npart[dir2]  = ngrid2;
    io::npart[dir3]  = ngrid3;
    
    io::nvars = 9;
    io::nvarg = 2;
    int ngrid1mod = ngrid1;
    int ngrid2mod = ngrid2;
    int ngrid3mod = ngrid3;
    // for grid, we include the edge
    if(coords[dir1] == ntiles[dir1] -1) ngrid1mod += 1;
    if(coords[dir2] == ntiles[dir2] -1) ngrid2mod += 1;
    if(coords[dir3] == ntiles[dir3] -1) ngrid3mod += 1;
    
    io::gridXout.allocate(io::nvarg,ngrid1mod);
    io::gridYout.allocate(io::nvarg,ngrid2mod);
    io::gridZout.allocate(io::nvarg,ngrid3mod);
    io::Fieldout.allocate(io::nvars,ngrid3,ngrid2,ngrid1);

    for(int i=0;i<=ngrid1mod;i++){
      io::gridXout(0,i) = G.x1b(is+i);
      io::gridXout(1,i) = G.x1a(is+i);
    }
    
    for(int j=0;j<=ngrid2mod;j++){
      io::gridYout(0,j) = G.x2b(js+j);
      io::gridYout(1,j) = G.x2a(js+j);
    }
    for(int k=0;k<=ngrid3mod;k++){
      io::gridZout(0,k) = G.x3b(ks+k);
      io::gridZout(1,k) = G.x3a(ks+k);
    }
    
    is_inited = true;
  };
  
  // ---- output text (unf%05d.dat) ----
  char fname_unf[256];
  std::snprintf(fname_unf, sizeof(fname_unf), "bindata/unf%05d.dat", index);
  FILE* fp_unf = std::fopen(fname_unf, "w");
  if (!fp_unf){
    std::fprintf(stderr, "open failed: %s : %s\n", fname_unf, std::strerror(errno));
  }
  
  std::fprintf(fp_unf, "# %.17g %.17g\n", time_sim,dt);
  std::fprintf(fp_unf, "# %d \n", ngrid1);
  std::fprintf(fp_unf, "# %d \n", ngrid2);
  std::fprintf(fp_unf, "# %d \n", ngrid3);
  std::fclose(fp_unf);
  
  // ---- output data (bin%05d.dat) ----

  for (int k=0;k<=ngrid3;k++)
    for (int j=0;j<=ngrid2;j++)
      for (int i=0;i<=ngrid1;i++){
	io::Fieldout(0,k,j,i) = P(nden,k+ks,j+js,i+is);
	io::Fieldout(1,k,j,i) = P(nve1,k+ks,j+js,i+is);
	io::Fieldout(2,k,j,i) = P(nve2,k+ks,j+js,i+is);
	io::Fieldout(3,k,j,i) = P(nve3,k+ks,j+js,i+is);
	io::Fieldout(4,k,j,i) = P(nbm1,k+ks,j+js,i+is);
	io::Fieldout(5,k,j,i) = P(nbm2,k+ks,j+js,i+is);
	io::Fieldout(6,k,j,i) = P(nbm3,k+ks,j+js,i+is);
	io::Fieldout(7,k,j,i) = P(nbps,k+ks,j+js,i+is);
	io::Fieldout(8,k,j,i) = P(npre,k+ks,j+js,i+is); 
  }
  io::MPIOutputBindary(index);
  index += 1;
}


int main() {
  
  using namespace resolution_mod;
  using namespace hydflux_mod;
  using namespace boundary_mod;
  using namespace mpi_config_mod;
  const bool NoOutput = false;
  static bool is_final = false;

  periodic[dir1] = 1;
  periodic[dir2] = 1;
  periodic[dir3] = 1;
  ntiles[dir1]   = 1;
  ntiles[dir2]   = 2;
  ntiles[dir3]   = 1;
  InitializeMPI();
  
  if(myid_w == 0) printf("setup grids and fields\n");
  
  AllocateHydroVariables(G,U,Fx,Fy,Fz,P);
 
  AllocateBoundaryVariables(Bs,Br);
  
  if (myid_w == 0) printf("grid size for x y z = %i %i %i\n",ngrid1*ntiles[dir1],ngrid2*ntiles[dir2],ngrid3*ntiles[dir3]);
  
  GenerateGrid(G);
  //GenerateProblem(G,P,U);
  GenerateProblem2(G,P,U);
  if (myid_w == 0) printf("entering main loop\n");
  int step = 0;
  auto time_beg = std::chrono::high_resolution_clock::now();

  for (step=0;step<stepmax;step++){
    ControlTimestep(G); 
    if (myid_w==0 && step%300 ==0 && !NoOutput) printf("step=%i time=%e dt=%e\n",step,time_sim,dt);
    //printf("step=%i time=%e dt=%e\n",step,time_sim,dt);
    SetBoundaryCondition(P,Bs,Br);
    EvaluateCh();
    GetNumericalFlux1(G,P,Fx);
    GetNumericalFlux2(G,P,Fy);
    GetNumericalFlux3(G,P,Fz);
    UpdateConservU(G,Fx,Fy,Fz,U);
    DampPsi(G,U);
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
  if (myid_w == 0) printf("exiting main loop time=%e, step=%i\n",time_sim,step);
  if (myid_w == 0) printf("sim time [s]: %e\n", elapsed.count());
  if (myid_w == 0) printf("time/count/cell : %e\n", elapsed.count()/(ngrid1*ngrid2*ngrid3)/stepmax);

  is_final = true;
  //Output(is_final);
  //Output1D(is_final);
  
  printf("program has been finished\n");
  return 0;
}
