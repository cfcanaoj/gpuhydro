/**
 * @file hydro.cpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <omp.h>

#include "resolution.hpp"
#include "hydro.hpp"

#include "mpi_config.hpp"

using namespace resolution_mod;

namespace hydflux_mod {
#pragma omp declare target
  GridArray<double> G;
  FieldArray<double> U; //! U(mconsv,ktot,jtot,itot)
  FieldArray<double> Fx,Fy,Fz;
  FieldArray<double> P; //! P(nprim,ktot,jtot,itot)
  double csiso;
  double chg;
#pragma omp end declare target

auto assoc = [&](void* host_ptr, size_t bytes, int dev) {
    void* dptr = omp_get_mapped_ptr(host_ptr, dev);
    if (!dptr) {
        std::fprintf(stderr, "mapped_ptr is NULL for %p\n", host_ptr);
    }
    int r = omp_target_associate_ptr(host_ptr, dptr, bytes, /*device_offset=*/0, dev);
    if (r != 0) {
        std::fprintf(stderr, "omp_target_associate_ptr failed (%d)\n", r);
    }
};

void AllocateHydroVariables(GridArray<double>& G,FieldArray<double>& U,FieldArray<double>& Fx,FieldArray<double>& Fy,FieldArray<double>& Fz,FieldArray<double>& P){
  
  int dev = omp_get_default_device();
  
  G.allocate(ktot,jtot,itot);
  
  for (int i=0; i<itot; i++) {
    G.x1a(i) = 0.0;
    G.x1b(i) = 0.0;
  }
  for (int j=0; j<jtot; j++) {
    G.x2a(j) = 0.0;
    G.x2b(j) = 0.0;
  }
  for (int k=0; k<ktot; k++) {
    G.x3a(k) = 0.0;
    G.x3b(k) = 0.0;
  }
  
#pragma omp target update to (G.n1, G.n2, G.n3)
#pragma omp target enter data map (alloc: G.x1a_data[0:G.n1],G.x1b_data[0:G.n1])
#pragma omp target enter data map (alloc: G.x2a_data[0:G.n2],G.x2b_data[0:G.n2])
#pragma omp target enter data map (alloc: G.x3a_data[0:G.n3],G.x3b_data[0:G.n3])
  /*
  assoc(G.x1a_data, sizeof(double)*G.n1,dev);
  assoc(G.x1b_data, sizeof(double)*G.n1,dev);
  assoc(G.x2a_data, sizeof(double)*G.n2,dev);
  assoc(G.x2b_data, sizeof(double)*G.n2,dev);
  assoc(G.x3a_data, sizeof(double)*G.n3,dev);
  assoc(G.x3b_data, sizeof(double)*G.n3,dev);
  */
  U.allocate(mconsv,ktot,jtot,itot);
  P.allocate(nprim ,ktot,jtot,itot);
      
  Fx.allocate(mconsv,ktot,jtot,itot);
  Fy.allocate(mconsv,ktot,jtot,itot);
  Fz.allocate(mconsv,ktot,jtot,itot);

#pragma omp target enter data map (alloc: U.data[0: U.size])
#pragma omp target enter data map (alloc:Fx.data[0:Fx.size])
#pragma omp target enter data map (alloc:Fy.data[0:Fy.size])
#pragma omp target enter data map (alloc:Fz.data[0:Fz.size])
#pragma omp target enter data map (alloc: P.data[0: P.size])
  /*
  assoc( U.data, sizeof(double)* U.size,dev);
  assoc(Fx.data, sizeof(double)*Fx.size,dev);
  assoc(Fy.data, sizeof(double)*Fy.size,dev);
  assoc(Fz.data, sizeof(double)*Fz.size,dev);
  assoc( P.data, sizeof(double)* P.size,dev);
  */
  for (int m=0; m<mconsv; m++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
  	   U(m,k,j,i) = 0.0;
	  Fx(m,k,j,i) = 0.0;
	  Fy(m,k,j,i) = 0.0;
	  Fz(m,k,j,i) = 0.0;
  }

  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  P(n,k,j,i) = 0.0;
  }
  
#pragma omp target update to ( U.data[0: U.size], U.n1, U.n2, U.n3, U.nv)
#pragma omp target update to (Fx.data[0:Fx.size],Fx.n1,Fx.n2,Fx.n3,Fx.nv)
#pragma omp target update to (Fy.data[0:Fy.size],Fy.n1,Fy.n2,Fy.n3,Fy.nv)
#pragma omp target update to (Fz.data[0:Fz.size],Fz.n1,Fz.n2,Fz.n3,Fz.nv)
#pragma omp target update to ( P.data[0: P.size], P.n1, P.n2, P.n3, P.nv)

}

void DeallocateHydroVariables(GridArray<double>& G,FieldArray<double>& U,FieldArray<double>& Fx,FieldArray<double>& Fy,FieldArray<double>& Fz,FieldArray<double>& P){

#pragma omp target exit data map (delete: G.x1a_data[0:G.n1],G.x1b_data[0:G.n1],G.x2a_data[0:G.n2],G.x2b_data[0:G.n2],G.x3a_data[0:G.n3],G.x3b_data[0:G.n3])
#pragma omp target exit data map (delete: U.data[0: U.size], U.n1, U.n2, U.n3, U.nv)
#pragma omp target exit data map (delete:Fx.data[0:Fx.size],Fx.n1,Fx.n2,Fx.n3,Fx.nv)
#pragma omp target exit data map (delete:Fy.data[0:Fy.size],Fy.n1,Fy.n2,Fy.n3,Fy.nv)
#pragma omp target exit data map (delete:Fz.data[0:Fz.size],Fz.n1,Fz.n2,Fz.n3,Fz.nv)
#pragma omp target exit data map (delete: P.data[0: P.size], P.n1, P.n2, P.n3, P.nv)

}

#pragma omp declare target
void vanLeer(const double& dsvp,const double& dsvm,double& dsv){
    if(dsvp * dsvm > 0.0e0){
      dsv = 2.0e0*dsvp *dsvm/(dsvp + dsvm);
    }else{
      dsv = 0.0e0;
    }
}
#pragma omp end declare target
#pragma omp declare target
void HLLD(const double (&leftst)[2*mconsv+madd],const double (&rigtst)[2*mconsv+madd],double (&nflux)[mconsv]){

//=====================================================================
//
// HLLD Scheme
//
// Purpose
// Calculation of Numerical Flux by HLLD method
//
// Reference
//
// Input
// Output
//=====================================================================

//----- U -----
// qql :: left state
// qqr :: right state
      double  rol,vxl,vyl,vzl,ptl,eel;
      double  ror,vxr,vyr,vzr,ptr,eer;
      double  rxl,ryl,rzl;
      double  rxr,ryr,rzr;
      double  bxs,byl,bzl;
      double      byr,bzr;
      double  ptst;

//----- U* ----
// qqlst ::  left state
// qqrst :: right state
      double  rolst,vxlst,vylst,vzlst,eelst;
      double  rorst,vxrst,vyrst,vzrst,eerst;
      double  rxlst,rylst,rzlst;
      double  rxrst,ryrst,rzrst;
      double        bylst,bzlst;
      double        byrst,bzrst;

//----- U** ----
// qqlst ::  left state
// qqrst :: right state
      double  vyldst,vzldst,eeldst;
      double  vyrdst,vzrdst,eerdst;
      double  ryldst,rzldst;
      double  ryrdst,rzrdst;
      double        byldst,bzldst;
      double        byrdst,bzrdst;

//----- flux ---
// fqql ::  left physical flux
// fqqr :: right physical flux
      double  frol,frxl,fryl,frzl,feel;
      double            fbyl,fbzl;
      double  fror,frxr,fryr,frzr,feer;
      double            fbyr,fbzr;

//----- wave speed ---
// sl ::  left-going fastest signal velocity
// sr :: right-going fastest signal velocity
// sm :: contact discontinuity velocity
// slst ::  left-going alfven velocity
// srst :: right-going alfven velocity
      double  sm,sl,sr,slst,srst;

// cfl :: left-state Fast wave velocity
// cfr :: right-sate Fast wave velocity
      double  cfl,cfr;

//--------------------
// temporary variables
      double  sdl,sdr,sdml,sdmr,isdml,isdmr,rosdl,rosdr;
      double  temp;
  
// no if
      double  sign1,maxs1,mins1;
      double  msl,msr,mslst,msrst,temp1,invsumro,sqrtror,sqrtrol,abbx;
      double  bxsq,temp_fst,eps,itf,vdbstl,vdbstr,signbx;

//----- Step 0. ----------------------------------------------------------|
      eps = 1.0e-30;
//---- Left state
        
        rol = leftst[mudn];
        eel = leftst[muet];
        rxl = leftst[muvu];
        ryl = leftst[muvv];
        rzl = leftst[muvw];
        vxl = leftst[muvu]/leftst[mudn];
        vyl = leftst[muvv]/leftst[mudn];
        vzl = leftst[muvw]/leftst[mudn];
        byl = leftst[mubv];
        bzl = leftst[mubw];
        ptl = leftst[mpre];

//---- Right state
        
        ror = rigtst[mudn];
        eer = rigtst[muet];
        rxr = rigtst[muvu];
        ryr = rigtst[muvv];
        rzr = rigtst[muvw];
        vxr = rigtst[muvu]/rigtst[mudn];
        vyr = rigtst[muvv]/rigtst[mudn];
        vzr = rigtst[muvw]/rigtst[mudn];
        byr = rigtst[mubv];
        bzr = rigtst[mubw];
        ptr = rigtst[mpre];

//----- Step 1. ----------------------------------------------------------|
// Compute wave left & right wave speed
//
        cfl = leftst[mcsp];
        cfr = rigtst[mcsp];

        sl = std::min(vxl,vxr)-std::max(cfl,cfr); // note sl is negative
        sr = std::max(vxl,vxr)+std::max(cfl,cfr);
//----- Step 2. ----------------------------------------------------------|
// compute L/R fluxs
//
// Left value
        frol = leftst[mfdn];
        feel = leftst[mfet];
        frxl = leftst[mfvu];
        fryl = leftst[mfvv];
        frzl = leftst[mfvw];
        fbyl = leftst[mfbv];
        fbzl = leftst[mfbw];

// Right value
        fror = rigtst[mfdn];
        feer = rigtst[mfet];
        frxr = rigtst[mfvu];
        fryr = rigtst[mfvv]; 
        frzr = rigtst[mfvw];
        fbyr = rigtst[mfbv];
        fbzr = rigtst[mfbw];


//----- Step 4. ----------------------------------------------------------|
// compute middle and alfven wave
//
        sdl = std::min(-1e-20,sl - vxl);
        sdr = std::max( 1e-20,sr - vxr);
        rosdl = rol*sdl;
        rosdr = ror*sdr;

        temp = 1.0e0/(rosdr - rosdl);
// Eq. 45
        sm = (rosdr*vxr - rosdl*vxl - ptr + ptl)*temp;
           
        sdml = std::min(-1e-20,sl - sm); isdml = 1.0e0/sdml;
        sdmr = std::max( 1e-20,sr - sm); isdmr = 1.0e0/sdmr;

//----- Step 5. ----------------------------------------------------------|
// compute intermediate states
//
// Eq. 49
        ptst = (rosdr*ptl-rosdl*ptr+rosdl*rosdr*(vxr-vxl))*temp;
		
//----- Step 5A. ----------------------------------------------------------|
// compute Ul*
//
           bxs = 0.5e0*(leftst[mubu]+rigtst[mubu]);
           bxsq = bxs*bxs;
           temp_fst = rosdl*sdml - bxsq;
           sign1 = std::copysign(1.0e0,std::abs(temp_fst)-eps);

           maxs1 = std::max(0.0e0,sign1);
           mins1 = std::min(0.0e0,sign1);

           itf = 1.0e0/(temp_fst+mins1);
           isdml = 1.0e0/sdml;

           temp = bxs*(sdl-sdml)*itf;
           rolst = maxs1*(rosdl*isdml) - mins1*rol;
           vxlst = maxs1*sm - mins1*vxl;
           rxlst = rolst*vxlst;
           
           vylst = maxs1*(vyl-byl*temp) - mins1*vyl;
           rylst = rolst*vylst;
           vzlst = maxs1*(vzl-bzl*temp) - mins1*vzl;
           rzlst = rolst*vzlst;
           
           temp = (rosdl*sdl-bxsq)*itf;
           bylst = maxs1*(byl*temp) - mins1*byl;
           bzlst = maxs1*(bzl*temp) - mins1*bzl;

           vdbstl = vxlst*bxs+vylst*bylst+vzlst*bzlst;
           eelst = maxs1*(sdl*eel - ptl*vxl + ptst*sm +      
               bxs*(vxl*bxs+vyl*byl+vzl*bzl-vdbstl))*isdml  
               - mins1*eel		;
           
//----- Step 5B. ----------------------------------------------------------|
// compute Ur*
//
           temp_fst = rosdr*sdmr - bxsq;
           sign1 = std::copysign(1.0e0,std::abs(temp_fst)-eps);
           maxs1 = std::max(0.0e0,sign1);
           mins1 = std::min(0.0e0,sign1);

           itf = 1.0e0/(temp_fst+mins1);
           isdmr = 1.0e0/sdmr;
           
           temp = bxs*(sdr-sdmr)*itf;
           rorst = maxs1*(rosdr*isdmr) - mins1*ror;
           vxrst = maxs1*sm - mins1*vxr;
           rxrst = rorst*vxrst;
           
           vyrst = maxs1*(vyr-byr*temp) - mins1*vyr;
           ryrst = rorst*vyrst;
           vzrst = maxs1*(vzr-bzr*temp) - mins1*vzr;
           rzrst = rorst*vzrst;
           
           temp = (rosdr*sdr-bxsq)*itf;
           byrst = maxs1*(byr*temp) - mins1*byr;
           bzrst = maxs1*(bzr*temp) - mins1*bzr;
				
           vdbstr = vxrst*bxs+vyrst*byrst+vzrst*bzrst;
           eerst = maxs1*((sdr*eer - ptr*vxr  + ptst*sm)*isdmr +   
               bxs*(vxr*bxs+vyr*byr+vzr*bzr-vdbstr)*isdmr)        
               - mins1*eer;
;
//----- Step 5C. ----------------------------------------------------------|;
// compute Ul** and Ur**;
//;
           sqrtrol = sqrt(rolst);
           sqrtror = sqrt(rorst);

           abbx = std::abs(bxs);
           signbx = std::copysign(1.0e0,bxs);
           sign1  = std::copysign(1.0e0,abbx-eps);

           maxs1 =  std::max(0e0,sign1);
           mins1 = -std::min(0e0,sign1);
           invsumro = maxs1/(sqrtrol + sqrtror);

           temp = invsumro*(sqrtrol*vylst + sqrtror*vyrst
               + signbx*(byrst-bylst));
           vyldst = vylst*mins1 + temp;
           vyrdst = vyrst*mins1 + temp;
           ryldst = rylst*mins1 + rolst * temp;
           ryrdst = ryrst*mins1 + rorst * temp;

           temp = invsumro*(sqrtrol*vzlst + sqrtror*vzrst  
               + signbx*(bzrst-bzlst));
           vzldst = vzlst*mins1 + temp;
           vzrdst = vzrst*mins1 + temp;
           rzldst = rzlst*mins1 + rolst * temp;
           rzrdst = rzrst*mins1 + rorst * temp;

           temp = invsumro*(sqrtrol*byrst + sqrtror*bylst
               + signbx*sqrtrol*sqrtror*(vyrst-vylst));
           byldst = bylst*mins1 + temp;
           byrdst = byrst*mins1 + temp;

           temp = invsumro*(sqrtrol*bzrst + sqrtror*bzlst
                + signbx*sqrtrol*sqrtror*(vzrst-vzlst));
           bzldst = bzlst*mins1 + temp;
           bzrdst = bzrst*mins1 + temp;
           
           temp = sm*bxs + vyldst*byldst + vzldst*bzldst;
           eeldst = eelst - sqrtrol*signbx*(vdbstl - temp)*maxs1;
           eerdst = eerst + sqrtror*signbx*(vdbstr - temp)*maxs1;
           
//----- Step 6. ----------------------------------------------------------|;
// compute flux;
           slst = (sm - abbx/sqrtrol)*maxs1;
           srst = (sm + abbx/sqrtror)*maxs1;

           sign1 =  std::copysign(1.0e0,sm);
           maxs1 =  std::max(0.0e0,sign1);
           mins1 = -std::min(0.0e0,sign1);

           msl = std::min(sl,0.0e0);
           msr = std::max(sr,0.0e0);
           mslst = std::min(slst,0.0e0);
           msrst = std::max(srst,0.0e0);

           temp = mslst-msl;
           temp1 = msrst-msr;

           nflux[mden] = (frol+(rolst-rol)*msl)*maxs1
                        +(fror+(rorst-ror)*msr)*mins1;
           nflux[meto] = (feel+(eelst-eel)*msl+(eeldst-eelst)*mslst)*maxs1
                        +(feer+(eerst-eer)*msr+(eerdst-eerst)*msrst)*mins1;
           nflux[mrvu] = (frxl+(rxlst-rxl)*msl)*maxs1
                        +(frxr+(rxrst-rxr)*msr)*mins1;
           nflux[mrvv] = (fryl+(rylst-ryl)*msl+(ryldst-rylst)*mslst)*maxs1
                        +(fryr+(ryrst-ryr)*msr+(ryrdst-ryrst)*msrst)*mins1;
           nflux[mrvw] = (frzl+(rzlst-rzl)*msl+(rzldst-rzlst)*mslst)*maxs1
                        +(frzr+(rzrst-rzr)*msr+(rzrdst-rzrst)*msrst)*mins1;
           nflux[mbmu] = 0e0;
           nflux[mbmv] = (fbyl+(bylst-byl)*msl+(byldst-bylst)*mslst)*maxs1
                        +(fbyr+(byrst-byr)*msr+(byrdst-byrst)*msrst)*mins1;
           nflux[mbmw] = (fbzl+(bzlst-bzl)*msl+(bzldst-bzlst)*mslst)*maxs1
                        +(fbzr+(bzrst-bzr)*msr+(bzrdst-bzrst)*msrst)*mins1;
}
#pragma omp end declare target

void GetNumericalFlux1(const GridArray<double>&G,const FieldArray<double>& P,FieldArray<double>& Fx){

#pragma omp target teams distribute parallel for collapse(2)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je; j++){
	for (int i=is; i<=ie+1; i++) {
	  // Pick up related cell
	  double Pleftc1[nprim];
	  double Pleftc2[nprim];
	  double Prigtc1[nprim];
	  double Prigtc2[nprim];  
	  for (int n=0; n<nprim; n++){
	    Pleftc1[n] = P(n,k,j,i-2);
	    Pleftc2[n] = P(n,k,j,i-1);
	    Prigtc1[n] = P(n,k,j,i  );
	    Prigtc2[n] = P(n,k,j,i+1);
	  }
	  //if(i==is) printf("Pleftc1: %e %e %e %e\n",Pleftc1[mden],Pleftc1[mrvu],Pleftc1[meto],Pleftc1[mubp]);
	  // Calculte Left state
	  double Plefte [nprim];
	  double dsvp,dsvm,dsv;
	  /* | Pleftc1   | Pleftc2 =>| Prigtc1   | Prigtc2   |  */
	  /*                     You are here                   */
	  for (int n=0; n<nprim; n++){
	    dsvp =  Prigtc1[n]- Pleftc2[n];
	    dsvm =              Pleftc2[n]- Pleftc1[n];
	    vanLeer(dsvp,dsvm,dsv);
	    Plefte[n] = Pleftc2[n] + 0.5e0*dsv;
	  }

	  //if(i==is)printf("Plefte: %e %e %e %e\n",Plefte[mden],Plefte[mrvu],Plefte[meto],Plefte[mubp]);
	  double Clefte [2*mconsv+madd];
	Clefte[mudn] = Plefte[nden]; // rho
	Clefte[muvu] = Plefte[nve1]*Plefte[nden]; // rho v_x
	Clefte[muvv] = Plefte[nve2]*Plefte[nden]; // rho v_y
	Clefte[muvw] = Plefte[nve3]*Plefte[nden]; // rho v_z
        Clefte[muet] = Plefte[nene]*Plefte[nden]  // e_i
	              +0.5e0*Plefte[nden]*(                  
                          +Plefte[nve1]*Plefte[nve1]                 
                          +Plefte[nve2]*Plefte[nve2]                 
                          +Plefte[nve3]*Plefte[nve3])  // + rho v^2/2
                      +0.5e0*             (                 
                          +Plefte[nbm1]*Plefte[nbm1]                 
                          +Plefte[nbm2]*Plefte[nbm2]                 
                          +Plefte[nbm3]*Plefte[nbm3]); // + B^2/2

	Clefte[mubu] = Plefte[nbm1]; // b_x
	Clefte[mubv] = Plefte[nbm2]; // b_y
	Clefte[mubw] = Plefte[nbm3]; // b_z
	Clefte[mubp] = Plefte[nbps]; // psi
	
        double  ptl = Plefte[npre] + ( Plefte[nbm1]*Plefte[nbm1]
                                      +Plefte[nbm2]*Plefte[nbm2]
				      +Plefte[nbm3]*Plefte[nbm3])/2.0e0;

	Clefte[mfdn] = Plefte[nden]*Plefte[nve1];
	Clefte[mfvu] = Plefte[nden]*Plefte[nve1]*Plefte[nve1] 
	                      + ptl-Plefte[nbm1]*Plefte[nbm1];
	Clefte[mfvv] = Plefte[nden]*Plefte[nve2]*Plefte[nve1]
	                           -Plefte[nbm2]*Plefte[nbm1];
        Clefte[mfvw] = Plefte[nden]*Plefte[nve3]*Plefte[nve1]
	                           -Plefte[nbm3]*Plefte[nbm1];
        Clefte[mfet] = (Clefte[muet]+ptl)*Plefte[nve1]
                           -( Plefte[nbm1]*Plefte[nve1]
                             +Plefte[nbm2]*Plefte[nve2]
			     +Plefte[nbm3]*Plefte[nve3])*Plefte[nbm1];

	Clefte[mfbu] =  0.0e0;
	Clefte[mfbv] =  Plefte[nbm2]*Plefte[nve1]
	               -Plefte[nve2]*Plefte[nbm1];
	Clefte[mfbw] =  Plefte[nbm3]*Plefte[nve1]
	               -Plefte[nve3]*Plefte[nbm1];
	Clefte[mfbp] = 0.0e0;  // psi
     
	double css = Plefte[ncsp]*Plefte[ncsp];
        double cts =  css // c_s^2*c_a^2;
	     + ( Plefte[nbm1]*Plefte[nbm1]  
                +Plefte[nbm2]*Plefte[nbm2]  
		+Plefte[nbm3]*Plefte[nbm3] )/Plefte[nden];

         Clefte[mcsp] = sqrt((cts +sqrt(cts*cts
                                  -4.0e0*css*Plefte[nbm1]*Plefte[nbm1] 
                                            /Plefte[nden])   
			      )/2.0e0);
         Clefte[mvel] = Plefte[nve1];//direction dependent
         Clefte[mpre] = ptl;
	 // Calculte Right state
	 
	  /* | Pleftc1   | Pleftc2 |<= Prigtc1   | Prigtc2   |  */
	  /*                     You are here                   */
	 double Prigte [nprim];
	  for (int n=0; n<nprim; n++){
	    dsvp =  Prigtc2[n]- Prigtc1[n];
	    dsvm =              Prigtc1[n]- Pleftc2[n];
	    vanLeer(dsvp,dsvm,dsv);
	    Prigte[n] = Prigtc1[n] - 0.5e0*dsv;
	  }
	  double Crigte [2*mconsv+madd];

	Crigte[mudn] = Prigte[nden]; // rho
	Crigte[muvu] = Prigte[nve1]*Prigte[nden]; // rho v_x
	Crigte[muvv] = Prigte[nve2]*Prigte[nden]; // rho v_y
	Crigte[muvw] = Prigte[nve3]*Prigte[nden]; // rho v_z
        Crigte[muet] = Prigte[nene]*Prigte[nden]  // e_i
	              +0.5e0*Prigte[nden]*(                  
                          +Prigte[nve1]*Prigte[nve1]                 
                          +Prigte[nve2]*Prigte[nve2]                 
                          +Prigte[nve3]*Prigte[nve3])  // + rho v^2/2
                      +0.5e0*             (                 
                          +Prigte[nbm1]*Prigte[nbm1]                 
                          +Prigte[nbm2]*Prigte[nbm2]                 
                          +Prigte[nbm3]*Prigte[nbm3]); // + B^2/2

	Crigte[mubu] = Prigte[nbm1]; // b_x
	Crigte[mubv] = Prigte[nbm2]; // b_y
	Crigte[mubw] = Prigte[nbm3]; // b_z
	Crigte[mubp] = Prigte[nbps]; // psi
	// total pressure
                 ptl = Prigte[npre] + ( Prigte[nbm1]*Prigte[nbm1]
                                       +Prigte[nbm2]*Prigte[nbm2]
				       +Prigte[nbm3]*Prigte[nbm3])/2.0e0;

	Crigte[mfdn] = Prigte[nden]*Prigte[nve1];
	Crigte[mfvu] = Prigte[nden]*Prigte[nve1]*Prigte[nve1] 
	                      + ptl-Prigte[nbm1]*Prigte[nbm1];
	Crigte[mfvv] = Prigte[nden]*Prigte[nve2]*Prigte[nve1]
	                           -Prigte[nbm2]*Prigte[nbm1];
        Crigte[mfvw] = Prigte[nden]*Prigte[nve3]*Prigte[nve1]
	                           -Prigte[nbm3]*Prigte[nbm1];
        Crigte[mfet] = (Crigte[muet]+ptl)*Prigte[nve1]
                           -( Prigte[nbm1]*Prigte[nve1]
                             +Prigte[nbm2]*Prigte[nve2]
			     +Prigte[nbm3]*Prigte[nve3])*Prigte[nbm1];

	Crigte[mfbu] =  0.0e0;
	Crigte[mfbv] =  Prigte[nbm2]*Prigte[nve1]
	               -Prigte[nve2]*Prigte[nbm1];
	Crigte[mfbw] =  Prigte[nbm3]*Prigte[nve1]
	               -Prigte[nve3]*Prigte[nbm1];
	Crigte[mfbp] = 0.0e0;  // psi
     
	       css = Prigte[ncsp]*Prigte[ncsp];
               cts =  css // c_s^2*c_a^2;
	     + ( Prigte[nbm1]*Prigte[nbm1]  
                +Prigte[nbm2]*Prigte[nbm2]  
		+Prigte[nbm3]*Prigte[nbm3] )/Prigte[nden];

         Crigte[mcsp] = sqrt((cts +sqrt(cts*cts
                                  -4.0e0*css*Prigte[nbm1]*Prigte[nbm1] 
                                            /Prigte[nden])   
			      )/2.0e0);
         Crigte[mvel] = Prigte[nve1]; //direction dependent
         Crigte[mpre] = ptl;
	 //Calculate flux
	 //if(i==is)printf("Cleft: %e %e %e %e\n",Clefte[mden],Clefte[mrvu],Clefte[meto],Clefte[mubp]);
	 double numflux [mconsv];
	 HLLD(Clefte,Crigte,numflux);
	 Fx(mden,k,j,i) = numflux[mden];
	 Fx(mrv1,k,j,i) = numflux[mrvu];
	 Fx(mrv2,k,j,i) = numflux[mrvv];
	 Fx(mrv3,k,j,i) = numflux[mrvw];
	 Fx(meto,k,j,i) = numflux[meto];
	 //Fx(mbm1,k,j,i) = numflux[mbmu];
	 Fx(mbm2,k,j,i) = numflux[mbmv];
	 Fx(mbm3,k,j,i) = numflux[mbmw];
	 
	 Fx(mbm1,k,j,i) = 0.5e0*    (Clefte[mubp]+Crigte[mubp])
	                 -0.5e0*chg*(Crigte[mubu]-Clefte[mubu]);
	 Fx(mbps,k,j,i) =(0.5e0*    (Clefte[mubu]+Crigte[mubu])
			  -0.5e0/chg*(Crigte[mubp]-Clefte[mubp]))*chg*chg;
	 
	}// i-loop
  }// j,k-loop
}
void GetNumericalFlux2(const GridArray<double>&G,const FieldArray<double>& P,FieldArray<double>& Fy){
  /* | Pleftc1   | Pleftc2 | Prigtc1   | Prigtc2   |              */
  /*                     You are here                             */
    
#pragma omp target teams distribute parallel for collapse(2)
    for (int k=ks; k<=ke; k++)
      for (int i=is; i<=ie; i++) {
	for (int j=js; j<=je+1; j++){
	  double Pleftc1[nprim];
	  double Pleftc2[nprim];
	  double Prigtc1[nprim];
	  double Prigtc2[nprim];	  
	  for (int n=0; n<nprim; n++){
	    Pleftc1[n] = P(n,k,j-2,i);
	    Pleftc2[n] = P(n,k,j-1,i);
	    Prigtc1[n] = P(n,k,j  ,i);
	    Prigtc2[n] = P(n,k,j+1,i);
	  }
	  double Plefte [nprim];
	  double dsvp,dsvm,dsv;
	  
	  /* | Pleftc1   | Pleftc2 =>| Prigtc1   | Prigtc2   |  */
	  /*                     You are here                   */
	  // Calculte Left state 
	  for (int n=0; n<nprim; n++){
	    dsvp =  Prigtc1[n]- Pleftc2[n];
	    dsvm =  Pleftc2[n]- Pleftc1[n];
	    vanLeer(dsvp,dsvm,dsv);
	    Plefte[n] = Pleftc2[n] + 0.5e0*dsv;
	  }
	  double Clefte [2*mconsv+madd];
	  //direction dependent w,u,v
	Clefte[mudn] = Plefte[nden]; // rho
	Clefte[muvw] = Plefte[nve1]*Plefte[nden]; // rho v_x
	Clefte[muvu] = Plefte[nve2]*Plefte[nden]; // rho v_y
	Clefte[muvv] = Plefte[nve3]*Plefte[nden]; // rho v_z
        Clefte[muet] = Plefte[nene]*Plefte[nden]  // e_i
	              +0.5e0*Plefte[nden]*(                  
                          +Plefte[nve1]*Plefte[nve1]                 
                          +Plefte[nve2]*Plefte[nve2]                 
                          +Plefte[nve3]*Plefte[nve3])  // + rho v^2/2
                      +0.5e0*             (                 
                          +Plefte[nbm1]*Plefte[nbm1]                 
                          +Plefte[nbm2]*Plefte[nbm2]                 
                          +Plefte[nbm3]*Plefte[nbm3]); // + B^2/2

	Clefte[mubw] = Plefte[nbm1]; // b_x
	Clefte[mubu] = Plefte[nbm2]; // b_y
	Clefte[mubv] = Plefte[nbm3]; // b_z
	Clefte[mubp] = Plefte[nbps]; // psi
	
        double  ptl = Plefte[npre] + ( Plefte[nbm1]*Plefte[nbm1]
                                      +Plefte[nbm2]*Plefte[nbm2]
				      +Plefte[nbm3]*Plefte[nbm3])/2.0e0;
	//direction dependent, nve2 or nbm2
	Clefte[mfdn] = Plefte[nden]*Plefte[nve2];
	Clefte[mfvw] = Plefte[nden]*Plefte[nve1]*Plefte[nve2] 
	                           -Plefte[nbm1]*Plefte[nbm2];
	Clefte[mfvu] = Plefte[nden]*Plefte[nve2]*Plefte[nve2]
      	                     + ptl -Plefte[nbm2]*Plefte[nbm2];// p diagnonal
        Clefte[mfvv] = Plefte[nden]*Plefte[nve3]*Plefte[nve2]
	                          - Plefte[nbm3]*Plefte[nbm2];
        Clefte[mfet] = (Clefte[muet]+ptl)*Plefte[nve2]
                           -( Plefte[nbm1]*Plefte[nve1]
                             +Plefte[nbm2]*Plefte[nve2]
			     +Plefte[nbm3]*Plefte[nve3])*Plefte[nbm2];

	Clefte[mfbw] =  Plefte[nbm1]*Plefte[nve2]
	              - Plefte[nve1]*Plefte[nbm2];
	Clefte[mfbu] =  0.0;
	Clefte[mfbv] =  Plefte[nbm3]*Plefte[nve2]
	              - Plefte[nve3]*Plefte[nbm2];
	Clefte[mfbp] = 0.0e0;  // psi
     
	double css = Plefte[ncsp]*Plefte[ncsp];
        double cts =  css // c_s^2*c_a^2;
	     + ( Plefte[nbm1]*Plefte[nbm1]  
                +Plefte[nbm2]*Plefte[nbm2]  
		+Plefte[nbm3]*Plefte[nbm3] )/Plefte[nden];

         Clefte[mcsp] = sqrt((cts +sqrt(cts*cts
                                  -4.0e0*css*Plefte[nbm2]*Plefte[nbm2]  //direction dependent
                                            /Plefte[nden])   
			      )/2.0e0);
         Clefte[mvel] = Plefte[nve2];//direction dependent
         Clefte[mpre] = ptl;
	 
	  /* | Pleftc1   | Pleftc2 |<= Prigtc1   | Prigtc2   |  */
	  /*                     You are here                   */
	 // Calculte Right state	
	  double Prigte [nprim];
	  for (int n=0; n<nprim; n++){
	    dsvp =  Prigtc2[n]- Prigtc1[n];
	    dsvm =  Prigtc1[n]- Pleftc2[n];
	    vanLeer(dsvp,dsvm,dsv);
	    Prigte[n] = Prigtc1[n] - 0.5e0*dsv;
	  }
	  //direction dependent w, u,v
	  double Crigte [2*mconsv+madd];
	Crigte[mudn] = Prigte[nden]; // rho
	Crigte[muvw] = Prigte[nve1]*Prigte[nden]; // rho v_x
	Crigte[muvu] = Prigte[nve2]*Prigte[nden]; // rho v_y
	Crigte[muvv] = Prigte[nve3]*Prigte[nden]; // rho v_z
        Crigte[muet] = Prigte[nene]*Prigte[nden]  // e_i
	              +0.5e0*Prigte[nden]*(                  
                          +Prigte[nve1]*Prigte[nve1]                 
                          +Prigte[nve2]*Prigte[nve2]                 
                          +Prigte[nve3]*Prigte[nve3])  // + rho v^2/2
                      +0.5e0*             (                 
                          +Prigte[nbm1]*Prigte[nbm1]                 
                          +Prigte[nbm2]*Prigte[nbm2]                 
                          +Prigte[nbm3]*Prigte[nbm3]); // + B^2/2

	Crigte[mubw] = Prigte[nbm1]; // b_x
	Crigte[mubu] = Prigte[nbm2]; // b_y
	Crigte[mubv] = Prigte[nbm3]; // b_z
	Crigte[mubp] = Prigte[nbps]; // psi
	// total pressure
                 ptl = Prigte[npre] + ( Prigte[nbm1]*Prigte[nbm1]
                                       +Prigte[nbm2]*Prigte[nbm2]
		          	       +Prigte[nbm3]*Prigte[nbm3])/2.0e0;

	//direction dependent, nve2 or nbm2
	Crigte[mfdn] = Prigte[nden]*Prigte[nve2];
	Crigte[mfvw] = Prigte[nden]*Prigte[nve1]*Prigte[nve2] 
	                           -Prigte[nbm1]*Prigte[nbm2];
	Crigte[mfvu] = Prigte[nden]*Prigte[nve2]*Prigte[nve2]
	                     + ptl -Prigte[nbm2]*Prigte[nbm2]; //diagnonal
        Crigte[mfvv] = Prigte[nden]*Prigte[nve3]*Prigte[nve2]
	                          - Prigte[nbm3]*Prigte[nbm2];
        Crigte[mfet] = (Crigte[muet]+ptl)*Prigte[nve2]
                           -( Prigte[nbm1]*Prigte[nve1]
                             +Prigte[nbm2]*Prigte[nve2]
			     +Prigte[nbm3]*Prigte[nve3])*Prigte[nbm2];

	Crigte[mfbw] =  Prigte[nbm1]*Prigte[nve2]
	               -Prigte[nve1]*Prigte[nbm2];
	Crigte[mfbu] =  0.0e0;
	Crigte[mfbv] =  Prigte[nbm3]*Prigte[nve2]
	               -Prigte[nve3]*Prigte[nbm2];
	Crigte[mfbp] = 0.0e0;  // psi
     
	       css = Prigte[ncsp]*Prigte[ncsp];
               cts =  css // c_s^2*c_a^2;
	     + ( Prigte[nbm1]*Prigte[nbm1]  
                +Prigte[nbm2]*Prigte[nbm2]  
		+Prigte[nbm3]*Prigte[nbm3] )/Prigte[nden];

         Crigte[mcsp] = sqrt((cts +sqrt(cts*cts
                                  -4.0e0*css*Prigte[nbm2]*Prigte[nbm2] 
                                            /Prigte[nden])   
			      )/2.0e0);
         Crigte[mvel] = Prigte[nve2];//direction dependent
         Crigte[mpre] = ptl;
	 //Calculate flux
	  double numflux [mconsv];  
	 HLLD(Clefte,Crigte,numflux);
	 Fy(mden,k,j,i) = numflux[mden];
	 Fy(mrv1,k,j,i) = numflux[mrvw];
	 Fy(mrv2,k,j,i) = numflux[mrvu];
	 Fy(mrv3,k,j,i) = numflux[mrvv];
	 Fy(meto,k,j,i) = numflux[meto];
	 Fy(mbm1,k,j,i) = numflux[mbmw];
	 //Fx(mbm2,k,j,i) = numflux[mbmu];
	 Fy(mbm3,k,j,i) = numflux[mbmv];
	 
	 Fy(mbm2,k,j,i) = 0.5e0*    (Clefte[mubp]+Crigte[mubp])
	                 -0.5e0*chg*(Crigte[mubu]-Clefte[mubu]);
	 Fy(mbps,k,j,i) =(0.5e0*    (Clefte[mubu]+Crigte[mubu])
	        	 -0.5e0/chg*(Crigte[mubp]-Clefte[mubp]))*chg*chg;
	 
	}// j-loop
  }// k,i-loop
}

void GetNumericalFlux3(const GridArray<double>&G,const FieldArray<double>& P,FieldArray<double>& Fz){
  
#pragma omp target teams distribute parallel for collapse(2)
  for (int j=js; j<=je; ++j)
    for (int i=is; i<=ie; ++i){
      for (int k=ks; k<=ke+1; ++k){
	double Pleftc1[nprim];
	double Pleftc2[nprim];
	double Prigtc1[nprim];
	double Prigtc2[nprim];
  
	  for (int n=0; n<nprim; n++){
	    Pleftc1[n] = P(n,k-2,j,i);
	    Pleftc2[n] = P(n,k-1,j,i);
	    Prigtc1[n] = P(n,k  ,j,i);
	    Prigtc2[n] = P(n,k+1,j,i);
	  }
	  
	  /* | Pleftc1   | Pleftc2 =>| Prigtc1   | Prigtc2   |  */
	  /*                     You are here                   */
	  // Calculte Left state
	  double Plefte [nprim];
	  double dsvp,dsvm,dsv;  
	  for (int n=0; n<nprim; n++){
	    dsvp =  Prigtc1[n]- Pleftc2[n];
	    dsvm =              Pleftc2[n]- Pleftc1[n];
	    vanLeer(dsvp,dsvm,dsv);
	    Plefte[n] = Pleftc2[n] + 0.5e0*dsv;
	  }
	  double Clefte [2*mconsv+madd];
	  //direction dependent v,w,u
	Clefte[mudn] = Plefte[nden]; // rho
	Clefte[muvv] = Plefte[nve1]*Plefte[nden]; // rho v_x
	Clefte[muvw] = Plefte[nve2]*Plefte[nden]; // rho v_y
	Clefte[muvu] = Plefte[nve3]*Plefte[nden]; // rho v_z
        Clefte[muet] = Plefte[nene]*Plefte[nden]  // e_i
	              +0.5e0*Plefte[nden]*(                  
                          +Plefte[nve1]*Plefte[nve1]                 
                          +Plefte[nve2]*Plefte[nve2]                 
                          +Plefte[nve3]*Plefte[nve3])  // + rho v^2/2
                      +0.5e0*             (                 
                          +Plefte[nbm1]*Plefte[nbm1]                 
                          +Plefte[nbm2]*Plefte[nbm2]                 
                          +Plefte[nbm3]*Plefte[nbm3]); // + B^2/2

	Clefte[mubv] = Plefte[nbm1]; // b_x
	Clefte[mubw] = Plefte[nbm2]; // b_y
	Clefte[mubu] = Plefte[nbm3]; // b_z
	Clefte[mubp] = Plefte[nbps]; // psi
	
        double  ptl = Plefte[npre] + ( Plefte[nbm1]*Plefte[nbm1]
                                      +Plefte[nbm2]*Plefte[nbm2]
				      +Plefte[nbm3]*Plefte[nbm3])/2.0e0;
	//direction dependent, nve3 or nbm3
	Clefte[mfdn] = Plefte[nden]*Plefte[nve3];
	Clefte[mfvv] = Plefte[nden]*Plefte[nve1]*Plefte[nve3] 
	                           -Plefte[nbm1]*Plefte[nbm3];
	Clefte[mfvw] = Plefte[nden]*Plefte[nve2]*Plefte[nve3]
      	                           -Plefte[nbm2]*Plefte[nbm3];
        Clefte[mfvu] = Plefte[nden]*Plefte[nve3]*Plefte[nve3]
	                     + ptl -Plefte[nbm3]*Plefte[nbm3]; // p diagnonal
        Clefte[mfet] = (Clefte[muet]+ptl)*Plefte[nve3]
                           -( Plefte[nbm1]*Plefte[nve1]
                             +Plefte[nbm2]*Plefte[nve2]
			     +Plefte[nbm3]*Plefte[nve3])*Plefte[nbm3];

	Clefte[mfbv] =  Plefte[nbm1]*Plefte[nve3]
	              - Plefte[nve1]*Plefte[nbm3];
	Clefte[mfbw] =  Plefte[nbm2]*Plefte[nve3]
	              - Plefte[nve2]*Plefte[nbm3];
	Clefte[mfbu] = 0.0e0;
	Clefte[mfbp] = 0.0e0;  // psi
     
	double css = Plefte[ncsp]*Plefte[ncsp];
        double cts =  css // c_s^2*c_a^2;
	     + ( Plefte[nbm1]*Plefte[nbm1]  
                +Plefte[nbm2]*Plefte[nbm2]  
		+Plefte[nbm3]*Plefte[nbm3] )/Plefte[nden];

         Clefte[mcsp] = sqrt((cts +sqrt(cts*cts
                                  -4.0e0*css*Plefte[nbm3]*Plefte[nbm3]  //direction dependent
                                            /Plefte[nden])   
			      )/2.0e0);
         Clefte[mvel] = Plefte[nve3];//direction
         Clefte[mpre] = ptl;

	  /* | Pleftc1   | Pleftc2 |<= Prigtc1   | Prigtc2   |  */
	  /*                     You are here                   */
	 // Calculte Right state
	  double Prigte [nprim];	  
	  for (int n=0; n<nprim; n++){
	    dsvp =  Prigtc2[n]- Prigtc1[n];
	    dsvm =              Prigtc1[n]- Pleftc2[n];
	    vanLeer(dsvp,dsvm,dsv);
	    Prigte[n] = Prigtc1[n] - 0.5e0*dsv;
	  }
	  //direction dependent w, u,v
	  double Crigte [2*mconsv+madd];
	Crigte[mudn] = Prigte[nden]; // rho
	Crigte[muvv] = Prigte[nve1]*Prigte[nden]; // rho v_x
	Crigte[muvw] = Prigte[nve2]*Prigte[nden]; // rho v_y
	Crigte[muvu] = Prigte[nve3]*Prigte[nden]; // rho v_z
        Crigte[muet] = Prigte[nene]*Prigte[nden]  // e_i
	              +0.5e0*Prigte[nden]*(                  
                          +Prigte[nve1]*Prigte[nve1]                 
                          +Prigte[nve2]*Prigte[nve2]                 
                          +Prigte[nve3]*Prigte[nve3])  // + rho v^2/2
                      +0.5e0*             (                 
                          +Prigte[nbm1]*Prigte[nbm1]                 
                          +Prigte[nbm2]*Prigte[nbm2]                 
                          +Prigte[nbm3]*Prigte[nbm3]); // + B^2/2

	Crigte[mubv] = Prigte[nbm1]; // b_x
	Crigte[mubw] = Prigte[nbm2]; // b_y
	Crigte[mubu] = Prigte[nbm3]; // b_z
	Crigte[mubp] = Prigte[nbps]; // psi
	// total pressure
                 ptl = Prigte[npre] + ( Prigte[nbm1]*Prigte[nbm1]
                                       +Prigte[nbm2]*Prigte[nbm2]
		          	       +Prigte[nbm3]*Prigte[nbm3])/2.0e0;

	//direction dependent, nve3 or nbm3
	Crigte[mfdn] = Prigte[nden]             *Prigte[nve3];
	Crigte[mfvv] = Prigte[nden]*Prigte[nve1]*Prigte[nve3] 
	                           -Prigte[nbm1]*Prigte[nbm3];
	Crigte[mfvw] = Prigte[nden]*Prigte[nve2]*Prigte[nve3]
	                           -Prigte[nbm2]*Prigte[nbm3]; 
        Crigte[mfvu] = Prigte[nden]*Prigte[nve3]*Prigte[nve3]
	                     + ptl -Prigte[nbm3]*Prigte[nbm3]; //diagnonal
        Crigte[mfet] = (Crigte[muet]       +ptl)*Prigte[nve3]
                           -( Prigte[nbm1]*Prigte[nve1]
                             +Prigte[nbm2]*Prigte[nve2]
			     +Prigte[nbm3]*Prigte[nve3])*Prigte[nbm3];

	Crigte[mfbv] =  Prigte[nbm1]*Prigte[nve3]
	               -Prigte[nve1]*Prigte[nbm3];
	Crigte[mfbw] =  Prigte[nbm2]*Prigte[nve3]
	               -Prigte[nve2]*Prigte[nbm3];
	Crigte[mfbu] = 0.0e0;
	Crigte[mfbp] = 0.0e0;  // psi
     
	       css = Prigte[ncsp]*Prigte[ncsp];
               cts =  css // c_s^2*c_a^2;
	     + ( Prigte[nbm1]*Prigte[nbm1]  
                +Prigte[nbm2]*Prigte[nbm2]  
		+Prigte[nbm3]*Prigte[nbm3] )/Prigte[nden];

         Crigte[mcsp] = sqrt((cts +sqrt(cts*cts
                                  -4.0e0*css*Prigte[nbm3]*Prigte[nbm3] 
                                            /Prigte[nden])   
			      )/2.0e0);
         Crigte[mvel] = Prigte[nve3];//direction dependent
         Crigte[mpre] = ptl;
	 //Calculate flux
	  double numflux [mconsv];
	 HLLD(Clefte,Crigte,numflux);
	 Fz(mden,k,j,i) = numflux[mden];
	 Fz(mrv1,k,j,i) = numflux[mrvv];
	 Fz(mrv2,k,j,i) = numflux[mrvw];
	 Fz(mrv3,k,j,i) = numflux[mrvu];
	 Fz(meto,k,j,i) = numflux[meto];
	 Fz(mbm1,k,j,i) = numflux[mbmv];
	 Fz(mbm2,k,j,i) = numflux[mbmw];
	 //Fz(mbm3,k,j,i) = numflux[mbm];
	 
	 Fz(mbm3,k,j,i) = 0.5e0*    (Clefte[mubp]+Crigte[mubp])
	                 -0.5e0*chg*(Crigte[mubu]-Clefte[mubu]);
	 Fz(mbps,k,j,i) =(0.5e0*    (Clefte[mubu]+Crigte[mubu])
	        	 -0.5e0/chg*(Crigte[mubp]-Clefte[mubp]))*chg*chg;
	 
	}// k-loop
  }// i,j-loop

}

void UpdateConservU(const GridArray<double>& G,const FieldArray<double>& Fx,const FieldArray<double>& Fy,const FieldArray<double>& Fz,FieldArray<double>& U){

  //printf("pre  U:%e %e %e %e\n", U(mden,ks,js,is), U(mrv3,ks,js,is), U(mbm3,ks,js,is), U(meto,ks,js,is));
  //printf("pre Fx:%e %e %e %e\n",Fx(mden,ks,js,is),Fx(mrv3,ks,js,is),Fx(mbm3,ks,js,is),Fx(meto,ks,js,is));
  //printf("pre Fy:%e %e %e %e\n",Fy(mden,ks,js,is),Fy(mrv3,ks,js,is),Fy(mbm3,ks,js,is),Fy(meto,ks,js,is));
  //printf("pre Fz:%e %e %e %e\n",Fz(mden,ks,js,is),Fz(mrv3,ks,js,is),Fz(mbm3,ks,js,is),Fz(meto,ks,js,is));
  
#pragma omp target teams distribute parallel for collapse(4)
  for (int m=0; m<mconsv; m++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je; j++)
	for (int i=is; i<=ie; i++) {
	  U(m,k,j,i) -= dt * ( (Fx(m,k,j,i+1) - Fx(m,k,j,i)) / (G.x1a(i+1)-G.x1a(i))
			      +(Fy(m,k,j+1,i) - Fy(m,k,j,i)) / (G.x2a(j+1)-G.x2a(j))
			      +(Fz(m,k+1,j,i) - Fz(m,k,j,i)) / (G.x3a(k+1)-G.x3a(k)) );
  }
  //printf("aft U:%e %e %e %e\n",U(mden,ks,js,is),U(mrv1,ks,js,is),U(mbm1,ks,js,is),U(meto,ks,js,is));
}


void UpdatePrimitvP(const FieldArray<double>& U,FieldArray<double>& P){

  //printf("U:et b1 b2 b3=%e %e %e %e\n",U(meto,ks,js,is),U(mbm1,ks,js,is),U(mbm2,ks,js,is),U(mbm3,ks,js,is));
#pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
         P(nden,k,j,i) = U(mden,k,j,i);
	 P(nve1,k,j,i) = U(mrv1,k,j,i)/U(mden,k,j,i);
         P(nve2,k,j,i) = U(mrv2,k,j,i)/U(mden,k,j,i);
         P(nve3,k,j,i) = U(mrv3,k,j,i)/U(mden,k,j,i);
	 double ekin = 0.5e0*( U(mrv1,k,j,i)*U(mrv1,k,j,i)
			      +U(mrv2,k,j,i)*U(mrv2,k,j,i)
			      +U(mrv3,k,j,i)*U(mrv3,k,j,i))/U(mden,k,j,i);
	 double emag = 0.5e0*( U(mbm1,k,j,i)*U(mbm1,k,j,i)
			      +U(mbm2,k,j,i)*U(mbm2,k,j,i)
			      +U(mbm3,k,j,i)*U(mbm3,k,j,i));
         P(nene,k,j,i) =  (U(meto,k,j,i)-ekin-emag)/U(mden,k,j,i);//specific internal energy
	 P(npre,k,j,i) =  U(mden,k,j,i) * csiso * csiso;
	 P(ncsp,k,j,i) =  csiso;
	 P(nbm1,k,j,i) =  U(mbm1,k,j,i);
	 P(nbm2,k,j,i) =  U(mbm2,k,j,i);
	 P(nbm3,k,j,i) =  U(mbm3,k,j,i);
	 P(nbps,k,j,i) =  U(mbps,k,j,i);
	 
      }
}

void ControlTimestep(const GridArray<double>& G){
  using namespace mpi_config_mod;
  const double huge = 1.0e90;
  double dtminl = huge;
#pragma omp target teams distribute parallel for reduction(min:dtminl) collapse(3)
  for (int k=ks; k<=ke; k++)
    for (int j=js; j<=je; j++)
      for (int i=is; i<=ie; i++) {
	double ctot = sqrt(     P(ncsp,k,j,i)*P(ncsp,k,j,i)
			   + (  P(nbm1,k,j,i)*P(nbm1,k,j,i)
			      + P(nbm2,k,j,i)*P(nbm2,k,j,i)
			      + P(nbm3,k,j,i)*P(nbm3,k,j,i) ) / P(nden,k,j,i)); 
	double dtminloc = std::min({ (G.x1a(i+1)-G.x1a(i))/(std::abs(P(nve1,k,j,i))+ctot)
				    ,(G.x2a(j+1)-G.x2a(j))/(std::abs(P(nve2,k,j,i))+ctot)
				    ,(G.x3a(k+1)-G.x3a(k))/(std::abs(P(nve3,k,j,i))+ctot)});
	dtminl = std::min(dtminloc,dtminl);
	/*  if(dtminloc < dtmin){//for debug
	  dtmin = dtminloc;
	  ip=i;
	  jp=j;
	  kp=k;
	  }*/
	
      }
  //  Here dtminl is in host
  double dtming;
  int myid_wg;
  MPIminfind(dtminl,myid_w,dtming,myid_wg);
  dt = 0.05e0*dtming;
#pragma omp target update to (dt)
}

void EvaluateCh(){
  using namespace mpi_config_mod;
  double chgloc = 0.0e0;
#pragma omp target teams distribute parallel for collapse(3) reduction(max:chgloc)
  for (int k=ks; k<=ke; k++)
    for (int j=js; j<=je; j++)
      for (int i=is; i<=ie; i++) {
	double css = P(ncsp,k,j,i)*P(ncsp,k,j,i);
	double ca1 = P(nbm1,k,j,i)*P(nbm1,k,j,i) / P(nden,k,j,i);
	double ca2 = P(nbm2,k,j,i)*P(nbm2,k,j,i) / P(nden,k,j,i);
	double ca3 = P(nbm3,k,j,i)*P(nbm3,k,j,i) / P(nden,k,j,i);
	double cts = css + ca1 + ca2 + ca3;
        double cm1 = sqrt((cts+sqrt(cts*cts-4.0e0*css*ca1))/2.0e0);
        double ch1 = (std::abs(P(nve1,k,j,i))+cm1);
        double cm2 = sqrt((cts+sqrt(cts*cts-4.0e0*css*ca2))/2.0e0);
        double ch2 = (std::abs(P(nve2,k,j,i))+cm2);
        double cm3 = sqrt((cts+sqrt(cts*cts-4.0e0*css*ca3))/2.0e0);
        double ch3 = (std::abs(P(nve3,k,j,i))+cm3);

        chgloc = std::max({chgloc,ch1,ch2,ch3});
      }
  // Here chgg is in host
  double chgg;
  int  myid_wg;
  MPImaxfind(chgloc,myid_w,chgg,myid_wg);
  chg = chgg;
#pragma omp target update to (chg)


}

void DampPsi(const GridArray<double>& G,FieldArray<double>& U){
  const double alphabp = 0.1e0;
#pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke; k++)
    for (int j=js; j<=je; j++)
      for (int i=is; i<=ie; i++) {
	double dhl = std::min({G.x1a(i+1)-G.x1a(i),G.x2a(j+1)-G.x2a(j),G.x3a(k+1)-G.x3a(k)});
	double taui = alphabp * chg/dhl;
	  U(mbps,k,j,i) = U(mbps,k,j,i) *(1.0e0-dt*taui);
      }
}


}; // end of namespace
