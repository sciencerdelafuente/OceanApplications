/**************************************************************************************/  
//  Reads a NetCDF file with oceanic data structured in an 3D Arakawa C-grid.
//  Computes particle trajectories with a fixed settling velocity.
//  It includes a map function that gives initial and final positions of particles
//  when traveling from one depth to a deeper one. 
//  In this case study, the input file is generated from ROMS in the Canary island basin.
//  (Mason, E. (2009). High-resolution modelling of the Canary Basin oceanic circulation.)
//  Related work: https://link.aps.org/doi/10.1103/PhysRevE.104.065111
/**************************************************************************************/  
//  Code written by Rebeca de la Fuente :  2020                                 
/**************************************************************************************/  


#include <iostream>
#include <netcdfcpp.h>
#include <string>
#include <netcdf.h>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string> 
#include <complex>
#include <cmath>
#include <math.h>
#include <bits/stdc++.h> 
#include <ctime>
#include <omp.h>


using namespace std;

#include "pvariables.h"

// g++ -o o s.cpp -lnetcdf_c++

double pi=acos(-1.0);
double rads=(pi/180.0);
double degrees=(180.0/pi);
double rearth=6371000.0; /* meters */
double secondsday=86400.0;

//Parameters
int initialt,timesim;
double secs_per_unit_time=172800; //2 days per unit time
double intstep_secs,intstep;
bool pintegration;
double sinking_velocity;

Variables point;

static const int NC_ERR = 2;
char ncfile[256];
char ncfileNetwork[256];
int xi_rho,eta_rho,xi_psi,eta_psi,xi_u,eta_v,s_rho,s_w,ntime,times,xi,eta,s;
double **lon_rho, **lat_rho, **h, **mask_rho,**lon_psi,**lat_psi;
double **lon_psi_meters, **lat_psi_meters,**mask_psi;

double ***free_surface,****u,****v,****w;
double ****depth_rho,****depth_w;
double ****depth_psi_w;

double ****up,****vp,****wp,****depthp;
int s_level=160;
double InitialPeriod;

int ReadingGridFlow();
int ReadingVelocityFlow();
int Scoordinate();

Hexaedron HX,HY,HZ;
double dxi,deta,ds;
double J[3][3];

int UVW(Variables &);
int Jacobian(Hexaedron &,Hexaedron &,Hexaedron &,double &,double &,double &);
double lon1,lat1,lon2,lat2;
double meters_u,meters_v,degrees_u,degrees_v;
double distance(double &,double &,double &,double &);
int conversion(double &,double &,double &,double &,double &,double &);
int conversion_uv();
double **relative_lon_psi_Q1_meters, **relative_lat_psi_Q1_meters, **relative_lon_psi_Q2_meters, **relative_lat_psi_Q2_meters, **relative_lon_psi_Q3_meters, **relative_lat_psi_Q3_meters;

int ParticleLocation(Variables &);
int ParticleLocation_auxiliar(Variables &);

bool InverseMapping_z(Variables &);
bool Boundaries(Variables &);
bool RK4(Variables &);
bool RK4_Backwards(Variables &);
int conversion_flow();
int VelocityFlow(Variables &);



int Grid();
int Down();
int FDLE_Depth();



int main(void)
{

   intstep_secs=20; 
   intstep=intstep_secs/secs_per_unit_time;


   InitialPeriod=390; 
   timesim=230; 
   sinking_velocity=-5/secondsday; //m/sec, negative sinking velocity = particle settles down from surface to bottom

   //-200 -> timesim=20 -> 16 GB
   //-30 -> timesim=40 -> 32 GB
   //-5 -> timesim=230 -> 128 GB

   ReadingGridFlow();
   ReadingVelocityFlow();
   Scoordinate();
   conversion_uv();
   conversion_flow();


   
   //las particulas se encuentran en la red con indices int 0- xi,eta,s,times


}





int FDLE_Depth()
{


      double GridResolution=0.5; // FDLE Resolution in xi,eta coordinates
      double InitialDepth=-100;
      double FinalDepth=-1200; 
      double time0=0; 

      double dx_depth=50;
      int number_files=-(FinalDepth-InitialDepth)/dx_depth;

      double time0days=time0*2;
      double Tmax=times-2;

    
      ofstream offile[number_files];
      string prefix,sufix;
      prefix="Sink390_VS5_";
      sufix=".txt";
      
      int kd=-InitialDepth+dx_depth;
      for(int i=0;i<number_files;i++)
      {
        stringstream ss;
        ss << prefix << kd << sufix;
        offile[i].open(ss.str().c_str());

      kd=kd+dx_depth;
      }




      double FinalDepthi;
      Variables particle;
      bool bc1,bc2,bc3,bc4;
      double anterior_xn,anterior_yn,anterior_lon,anterior_lat,anterior_depth,anterior_tau;
      double posterior_xn,posterior_yn,posterior_lon,posterior_lat,posterior_depth,posterior_tau;
      double lon_fixed,lat_fixed,time_fixed,depth_fixed,xi_fixed,eta_fixed;

      double xi0,eta0,s0,lon0,lat0,depth0;
      double lambda;

      double ktime;
      int particles_boundaries;
      int ncross;


      int contador=0;
      for(double part_xi=2;part_xi<xi-2;part_xi=part_xi+GridResolution)
      {
      for(double part_eta=2;part_eta<eta-2;part_eta=part_eta+GridResolution)
      {





          particle.t=time0;
          particle.xi=part_xi;
          particle.eta=part_eta;
          particle.depth=InitialDepth;
          bc1=InverseMapping_z(particle);
          bc2=Boundaries(particle);
          particle.domain=bc1*bc2;



          if(particle.domain==1)
          {
             ParticleLocation(particle);        
          }

          xi0=particle.xi;
          eta0=particle.eta;
          s0=particle.s;
          lon0=particle.lon;
          lat0=particle.lat;
          depth0=particle.depth;


          FinalDepthi=InitialDepth-dx_depth;
          ncross=0;
          ktime=time0;
          do{


                  ParticleLocation(particle); 
                  anterior_xn=particle.xi;
                  anterior_yn=particle.eta;
                  anterior_lon=particle.lon;
                  anterior_lat=particle.lat;
                  anterior_depth=particle.depth;
                  anterior_tau=particle.time_days;


                  bc1=RK4(particle);
                  bc2=Boundaries(particle);
                  particle.domain=particle.domain*bc1*bc2;

                  if(particle.domain==1)
                  {
                  ParticleLocation(particle); 
                  posterior_xn=particle.xi;
                  posterior_yn=particle.eta;
                  posterior_lon=particle.lon;
                  posterior_lat=particle.lat;
                  posterior_depth=particle.depth;
                  posterior_tau=particle.time_days;



                     if(anterior_depth>=FinalDepthi && posterior_depth<=FinalDepthi) //bajando
                     {
                        lambda=(FinalDepthi-anterior_depth)/(posterior_depth-anterior_depth);
                        lon_fixed=anterior_lon+lambda*(posterior_lon-anterior_lon);
                        lat_fixed=anterior_lat+lambda*(posterior_lat-anterior_lat);
                        time_fixed=anterior_tau+lambda*(posterior_tau-anterior_tau);
                        depth_fixed=anterior_depth+lambda*(posterior_depth-anterior_depth);
                        xi_fixed=anterior_xn+lambda*(posterior_xn-anterior_xn);
                        eta_fixed=anterior_yn+lambda*(posterior_yn-anterior_yn);

               offile[ncross]<<contador<<" "<<xi0<<" "<<eta0<<" "<<depth0<<" "<<lon0<<" "<<lat0<<" "<<xi_fixed<<" "<<eta_fixed<<" "<<depth_fixed<<" "<<lon_fixed<<" "<<lat_fixed<<" "<<time_fixed<<" "<<(time_fixed-time0days)<<endl;

                        ncross++;
                        FinalDepthi=FinalDepthi-dx_depth;
                      

                     }
                }


          ktime=ktime+intstep;
          }while(ktime<Tmax && particle.domain==1 && ncross<number_files);




         
         

      contador++;

      }
      }


}































int Down()
{


      double GridResolution=0.5; //FDLE Resolution in xi,eta coordinates
      double InitialDepth=-100;
    //  double FinalDepth0=InitialDepth-50;
    //  double FinalDepth1=InitialDepth-100;
      double FinalDepth2=-300;
      double time0=0; 



      double time0days=time0*2;
      double Tmax=times-2;
      Variables particle;
      bool bc1,bc2,bc3,bc4;
      int cross,cross0,cross1,cross2;

    
      ofstream offile0,offile1,offile2;
     // string filename0 = "05Down_z100_z150_vs5_";
     // string filename1 = "05Down_z100_z200_vs5";
      string filename2 = "05Down_z100_z300_vs200_St525";

     // offile0.open(filename0.c_str());
     // offile1.open(filename1.c_str());
      offile2.open(filename2.c_str());




      double anterior_xn,anterior_yn,anterior_lon,anterior_lat,anterior_depth,anterior_tau;
      double posterior_xn,posterior_yn,posterior_lon,posterior_lat,posterior_depth,posterior_tau;
      double lon_fixed,lat_fixed,time_fixed,depth_fixed,xi_fixed,eta_fixed;

      double xi0,eta0,s0,lon0,lat0,depth0;
      double lambda;

      double ktime;
      int particles_boundaries;

      int contador=0;
      for(double part_xi=2;part_xi<xi-2;part_xi=part_xi+GridResolution)
      {
      for(double part_eta=2;part_eta<eta-2;part_eta=part_eta+GridResolution)
      {





          particle.t=time0;
          particle.xi=part_xi;
          particle.eta=part_eta;
          particle.depth=InitialDepth;
          bc1=InverseMapping_z(particle);
          bc2=Boundaries(particle);
          particle.domain=bc1*bc2;



          if(particle.domain==1)
          {
             ParticleLocation(particle);        
          }

          xi0=particle.xi;
          eta0=particle.eta;
          s0=particle.s;
          lon0=particle.lon;
          lat0=particle.lat;
          depth0=particle.depth;


          cross=0;
          cross0=0; 
          cross1=0; 
          cross2=0; 
          ktime=time0;
          do{


                  ParticleLocation(particle); 
                  anterior_xn=particle.xi;
                  anterior_yn=particle.eta;
                  anterior_lon=particle.lon;
                  anterior_lat=particle.lat;
                  anterior_depth=particle.depth;
                  anterior_tau=particle.time_days;


                  bc1=RK4(particle);
                  bc2=Boundaries(particle);
                  particle.domain=particle.domain*bc1*bc2;

                  if(particle.domain==1)
                  {
                  ParticleLocation(particle); 
                  posterior_xn=particle.xi;
                  posterior_yn=particle.eta;
                  posterior_lon=particle.lon;
                  posterior_lat=particle.lat;
                  posterior_depth=particle.depth;
                  posterior_tau=particle.time_days;



                     if(anterior_depth>=FinalDepth2 && posterior_depth<=FinalDepth2 && cross2==0) 
                     {
                        lambda=(FinalDepth2-anterior_depth)/(posterior_depth-anterior_depth);
                        lon_fixed=anterior_lon+lambda*(posterior_lon-anterior_lon);
                        lat_fixed=anterior_lat+lambda*(posterior_lat-anterior_lat);
                        time_fixed=anterior_tau+lambda*(posterior_tau-anterior_tau);
                        depth_fixed=anterior_depth+lambda*(posterior_depth-anterior_depth);
                        xi_fixed=anterior_xn+lambda*(posterior_xn-anterior_xn);
                        eta_fixed=anterior_yn+lambda*(posterior_yn-anterior_yn);

               offile2<<contador<<" "<<xi0<<" "<<eta0<<" "<<depth0<<" "<<lon0<<" "<<lat0<<" "<<xi_fixed<<" "<<eta_fixed<<" "<<depth_fixed<<" "<<lon_fixed<<" "<<lat_fixed<<" "<<time_fixed<<" "<<(time_fixed-time0days)<<endl;

                        cross2=1;

                     }
                }

                cross=cross2;

          ktime=ktime+intstep;
          }while(ktime<Tmax && particle.domain==1 && cross==0);




         
         

      contador++;

      }
      }


}











































int Grid()
{


      double GridResolution=0.5; //FDLE Resolution in xi,eta coordinates
      double InitialDepth=-100;
      double FinalDepth0=InitialDepth-10;
      double FinalDepth1=InitialDepth-50;
      double FinalDepth2=InitialDepth-100;
      double time0=0; 



      double time0days=time0*2;
      double Tmax=times-2;
      Variables particle;
      bool bc1,bc2,bc3,bc4;
      int cross,cross0,cross1,cross2;

    


      double anterior_xn,anterior_yn,anterior_lon,anterior_lat,anterior_depth,anterior_tau;
      double posterior_xn,posterior_yn,posterior_lon,posterior_lat,posterior_depth,posterior_tau;
      double lon_fixed,lat_fixed,time_fixed,depth_fixed,xi_fixed,eta_fixed;

      double xi0,eta0,s0,lon0,lat0,depth0;
      double lambda;

      double ktime;
      int particles_boundaries;

      int contador=0;
      int intxi=0;
      for(double part_xi=2;part_xi<xi-2;part_xi=part_xi+GridResolution)
      {
      int inteta=0;
      for(double part_eta=2;part_eta<eta-2;part_eta=part_eta+GridResolution)
      {



         cout<<contador<<" "<<intxi<<" "<<inteta<<" "<<part_xi<<" "<<part_eta<<endl;
         

      contador++;

      inteta++;
      }
      intxi++;
      }


}















bool RK4_Backwards(Variables &point)  // RK4(xi,eta,t,s)
{

     bool inside=false;
     double K1x,K2x,K3x,K4x;
     double K1y,K2y,K3y,K4y;
     double K1z,K2z,K3z,K4z;
     Variables point1,point2,point3,point4;
     bool bc1,bc2,bc3,bc4;

     //1
     point1=point;

     bc1=Boundaries(point1);
     if(bc1==1)
     {
     VelocityFlow(point1);
     K1x=point1.u;
     K1y=point1.v;
     K1z=point1.w;
     

     //2
     point2.xi=point1.xi-0.5*K1x*intstep;
     point2.eta=point1.eta-0.5*K1y*intstep;
     point2.s=point1.s-0.5*K1z*intstep;
     point2.t=point1.t-0.5*intstep;

     bc2=Boundaries(point2);
     if(bc2==1)
     {
     VelocityFlow(point2);
     K2x=point2.u;
     K2y=point2.v;
     K2z=point2.w;

     //3
     point3.xi=point1.xi-0.5*K2x*intstep;
     point3.eta=point1.eta-0.5*K2y*intstep;
     point3.s=point1.s-0.5*K2z*intstep;
     point3.t=point1.t-0.5*intstep;

     bc3=Boundaries(point3);
     if(bc3==1)
     {
     VelocityFlow(point3);
     K3x=point3.u;
     K3y=point3.v;
     K3z=point3.w;

     //4
     point4.xi=point1.xi-K3x*intstep;
     point4.eta=point1.eta-K3y*intstep;
     point4.s=point1.s-K3z*intstep;
     point4.t=point1.t-intstep;

     bc4=Boundaries(point4);
     if(bc4==1)
     {
     VelocityFlow(point4);
     K4x=point4.u;
     K4y=point4.v;
     K4z=point4.w;


         point.xi=point.xi-intstep*((K1x+2*K2x+2*K3x+K4x)/6);
         point.eta=point.eta-intstep*((K1y+2*K2y+2*K3y+K4y)/6);
         point.s=point.s-intstep*((K1z+2*K2z+2*K3z+K4z)/6);
         point.t=point.t-intstep;

         inside=true;


     }
     }
     }
     }


     return inside;
}







bool RK4(Variables &point)  // RK4(xi,eta,t,s)
{

     bool inside=false;
     double K1x,K2x,K3x,K4x;
     double K1y,K2y,K3y,K4y;
     double K1z,K2z,K3z,K4z;
     Variables point1,point2,point3,point4;
     bool bc1,bc2,bc3,bc4;

     //1
     point1=point;

     bc1=Boundaries(point1);
     if(bc1==1)
     {
     VelocityFlow(point1);
     K1x=point1.u;
     K1y=point1.v;
     K1z=point1.w;
     

     //2
     point2.xi=point1.xi+0.5*K1x*intstep;
     point2.eta=point1.eta+0.5*K1y*intstep;
     point2.s=point1.s+0.5*K1z*intstep;
     point2.t=point1.t+0.5*intstep;

     bc2=Boundaries(point2);
     if(bc2==1)
     {
     VelocityFlow(point2);
     K2x=point2.u;
     K2y=point2.v;
     K2z=point2.w;

     //3
     point3.xi=point1.xi+0.5*K2x*intstep;
     point3.eta=point1.eta+0.5*K2y*intstep;
     point3.s=point1.s+0.5*K2z*intstep;
     point3.t=point1.t+0.5*intstep;

     bc3=Boundaries(point3);
     if(bc3==1)
     {
     VelocityFlow(point3);
     K3x=point3.u;
     K3y=point3.v;
     K3z=point3.w;

     //4
     point4.xi=point1.xi+K3x*intstep;
     point4.eta=point1.eta+K3y*intstep;
     point4.s=point1.s+K3z*intstep;
     point4.t=point1.t+intstep;

     bc4=Boundaries(point4);
     if(bc4==1)
     {
     VelocityFlow(point4);
     K4x=point4.u;
     K4y=point4.v;
     K4z=point4.w;


         point.xi=point.xi+intstep*((K1x+2*K2x+2*K3x+K4x)/6);
         point.eta=point.eta+intstep*((K1y+2*K2y+2*K3y+K4y)/6);
         point.s=point.s+intstep*((K1z+2*K2z+2*K3z+K4z)/6);
         point.t=point.t+intstep;

         inside=true;


     }
     }
     }
     }


     return inside;
}


















bool InverseMapping_z(Variables &point)  //xi,eta,t,depth -> s
{

    Variables point1,point2;
    double fx;
    bool inside=false;


    for(int i=0;i<s-1;i++)
    {

       point1.xi=point.xi;
       point1.eta=point.eta;
       point1.s=i;
       point1.t=point.t;

       ParticleLocation(point1);

       point2.xi=point.xi;
       point2.eta=point.eta;
       point2.s=i+1;
       point2.t=point.t;
       ParticleLocation(point2);

       if(point.depth>=point1.depth && point.depth<point2.depth)
       {
           fx=(i)+((point.depth-point1.depth)/(point2.depth-point1.depth));
           inside=true;
       }
    }

    point.s=fx;
    return inside;

}





int VelocityFlow(Variables &point) //xi,eta,s,t -> lon,lat,depth,time_days
{

   double Q11,Q12,Q21,Q22;
   double dt;

  // int psint=point.s;

   double ps=((point.s)*(s_level-1))/(s_w-1);

   //psi points
   int xi_i=point.xi;
   int eta_i=point.eta;
   int s_i=ps;
   int t_i=point.t;

   double x=point.xi-xi_i;
   double y=point.eta-eta_i;
   double sg=ps-s_i;
   double tg=point.t-t_i;

   dxi=x;
   deta=y;
   ds=sg;
   dt=tg;

   double u_s0_t0,u_s0_t1,u_s1_t0,u_s1_t1;
   double u_t0,u_t1;
   double v_s0_t0,v_s0_t1,v_s1_t0,v_s1_t1;
   double v_t0,v_t1;
   double w_s0_t0,w_s0_t1,w_s1_t0,w_s1_t1;
   double w_t0,w_t1;


  //U


  Q11=up[t_i][s_i][eta_i][xi_i];
  Q21=up[t_i][s_i][eta_i][xi_i+1];
  Q12=up[t_i][s_i][eta_i+1][xi_i];
  Q22=up[t_i][s_i][eta_i+1][xi_i+1];
  u_s0_t0=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=up[t_i][s_i+1][eta_i][xi_i];
  Q21=up[t_i][s_i+1][eta_i][xi_i+1];
  Q12=up[t_i][s_i+1][eta_i+1][xi_i];
  Q22=up[t_i][s_i+1][eta_i+1][xi_i+1];
  u_s1_t0=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=up[t_i+1][s_i][eta_i][xi_i];
  Q21=up[t_i+1][s_i][eta_i][xi_i+1];
  Q12=up[t_i+1][s_i][eta_i+1][xi_i];
  Q22=up[t_i+1][s_i][eta_i+1][xi_i+1];
  u_s0_t1=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=up[t_i+1][s_i+1][eta_i][xi_i];
  Q21=up[t_i+1][s_i+1][eta_i][xi_i+1];
  Q12=up[t_i+1][s_i+1][eta_i+1][xi_i];
  Q22=up[t_i+1][s_i+1][eta_i+1][xi_i+1];
  u_s1_t1=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  u_t0=sg*u_s1_t0+(1-sg)*u_s0_t0;
  u_t1=sg*u_s1_t1+(1-sg)*u_s0_t1;
  point.u=u_t1*tg+(1-tg)*u_t0;


  //V


  Q11=vp[t_i][s_i][eta_i][xi_i];
  Q21=vp[t_i][s_i][eta_i][xi_i+1];
  Q12=vp[t_i][s_i][eta_i+1][xi_i];
  Q22=vp[t_i][s_i][eta_i+1][xi_i+1];
  v_s0_t0=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=vp[t_i][s_i+1][eta_i][xi_i];
  Q21=vp[t_i][s_i+1][eta_i][xi_i+1];
  Q12=vp[t_i][s_i+1][eta_i+1][xi_i];
  Q22=vp[t_i][s_i+1][eta_i+1][xi_i+1];
  v_s1_t0=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=vp[t_i+1][s_i][eta_i][xi_i];
  Q21=vp[t_i+1][s_i][eta_i][xi_i+1];
  Q12=vp[t_i+1][s_i][eta_i+1][xi_i];
  Q22=vp[t_i+1][s_i][eta_i+1][xi_i+1];
  v_s0_t1=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=vp[t_i+1][s_i+1][eta_i][xi_i];
  Q21=vp[t_i+1][s_i+1][eta_i][xi_i+1];
  Q12=vp[t_i+1][s_i+1][eta_i+1][xi_i];
  Q22=vp[t_i+1][s_i+1][eta_i+1][xi_i+1];
  v_s1_t1=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  v_t0=sg*v_s1_t0+(1-sg)*v_s0_t0;
  v_t1=sg*v_s1_t1+(1-sg)*v_s0_t1;
  point.v=v_t1*tg+(1-tg)*v_t0;

  //W

  Q11=wp[t_i][s_i][eta_i][xi_i];
  Q21=wp[t_i][s_i][eta_i][xi_i+1];
  Q12=wp[t_i][s_i][eta_i+1][xi_i];
  Q22=wp[t_i][s_i][eta_i+1][xi_i+1];
  w_s0_t0=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=wp[t_i][s_i+1][eta_i][xi_i];
  Q21=wp[t_i][s_i+1][eta_i][xi_i+1];
  Q12=wp[t_i][s_i+1][eta_i+1][xi_i];
  Q22=wp[t_i][s_i+1][eta_i+1][xi_i+1];
  w_s1_t0=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=wp[t_i+1][s_i][eta_i][xi_i];
  Q21=wp[t_i+1][s_i][eta_i][xi_i+1];
  Q12=wp[t_i+1][s_i][eta_i+1][xi_i];
  Q22=wp[t_i+1][s_i][eta_i+1][xi_i+1];
  w_s0_t1=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=wp[t_i+1][s_i+1][eta_i][xi_i];
  Q21=wp[t_i+1][s_i+1][eta_i][xi_i+1];
  Q12=wp[t_i+1][s_i+1][eta_i+1][xi_i];
  Q22=wp[t_i+1][s_i+1][eta_i+1][xi_i+1];
  w_s1_t1=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  w_t0=sg*w_s1_t0+(1-sg)*w_s0_t0;
  w_t1=sg*w_s1_t1+(1-sg)*w_s0_t1;
  point.w=w_t1*tg+(1-tg)*w_t0;


}







int ParticleLocation(Variables &point) //xi,eta,s,t -> lon,lat,depth,time_days
{


   double Q11,Q12,Q21,Q22;
   double depth_s0_t0,depth_s0_t1,depth_s1_t0,depth_s1_t1;
   double depth_t0,depth_t1;
   double dt;

   //psi points
   int xi_i=point.xi;
   int eta_i=point.eta;
   int s_i=point.s;
   int t_i=point.t;

   double x=point.xi-xi_i;
   double y=point.eta-eta_i;
   double sg=point.s-s_i;
   double tg=point.t-t_i;

   dxi=x;
   deta=y;
   ds=sg;
   dt=tg;

  //time

  point.time_days=(secs_per_unit_time*point.t)/secondsday;

  //LON,LAT

  Q11=lon_psi[eta_i][xi_i];
  Q21=lon_psi[eta_i][xi_i+1];
  Q12=lon_psi[eta_i+1][xi_i];
  Q22=lon_psi[eta_i+1][xi_i+1];
  point.lon=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=lat_psi[eta_i][xi_i];
  Q21=lat_psi[eta_i][xi_i+1];
  Q12=lat_psi[eta_i+1][xi_i];
  Q22=lat_psi[eta_i+1][xi_i+1];
  point.lat=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));


  //Depth

  if(point.s<s_w-1)
  {


  Q11=depth_psi_w[t_i][s_i][eta_i][xi_i];
  Q21=depth_psi_w[t_i][s_i][eta_i][xi_i+1];
  Q12=depth_psi_w[t_i][s_i][eta_i+1][xi_i];
  Q22=depth_psi_w[t_i][s_i][eta_i+1][xi_i+1];
  depth_s0_t0=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=depth_psi_w[t_i][s_i+1][eta_i][xi_i];
  Q21=depth_psi_w[t_i][s_i+1][eta_i][xi_i+1];
  Q12=depth_psi_w[t_i][s_i+1][eta_i+1][xi_i];
  Q22=depth_psi_w[t_i][s_i+1][eta_i+1][xi_i+1];
  depth_s1_t0=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=depth_psi_w[t_i+1][s_i][eta_i][xi_i];
  Q21=depth_psi_w[t_i+1][s_i][eta_i][xi_i+1];
  Q12=depth_psi_w[t_i+1][s_i][eta_i+1][xi_i];
  Q22=depth_psi_w[t_i+1][s_i][eta_i+1][xi_i+1];
  depth_s0_t1=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  Q11=depth_psi_w[t_i+1][s_i+1][eta_i][xi_i];
  Q21=depth_psi_w[t_i+1][s_i+1][eta_i][xi_i+1];
  Q12=depth_psi_w[t_i+1][s_i+1][eta_i+1][xi_i];
  Q22=depth_psi_w[t_i+1][s_i+1][eta_i+1][xi_i+1];
  depth_s1_t1=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));

  depth_t0=sg*depth_s1_t0+(1-sg)*depth_s0_t0;
  depth_t1=sg*depth_s1_t1+(1-sg)*depth_s0_t1;
  point.depth=depth_t1*tg+(1-tg)*depth_t0;
 }

 if(point.s==s_w-1)
 {

  Q11=depth_psi_w[t_i][s-1][eta_i][xi_i];
  Q21=depth_psi_w[t_i][s-1][eta_i][xi_i+1];
  Q12=depth_psi_w[t_i][s-1][eta_i+1][xi_i];
  Q22=depth_psi_w[t_i][s-1][eta_i+1][xi_i+1];
  depth_s0_t0=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));


  Q11=depth_psi_w[t_i+1][s-1][eta_i][xi_i];
  Q21=depth_psi_w[t_i+1][s-1][eta_i][xi_i+1];
  Q12=depth_psi_w[t_i+1][s-1][eta_i+1][xi_i];
  Q22=depth_psi_w[t_i+1][s-1][eta_i+1][xi_i+1];
  depth_s0_t1=(Q11*(1-x)*(1-y))+(Q21*(x)*(1-y))+(Q12*(1-x)*(y))+(Q22*(x)*(y));


  point.depth=depth_s0_t1*tg+(1-tg)*depth_s0_t0;

 }

}








int conversion_flow()
{
  up = new double ***[times];  
  vp = new double ***[times]; 
  wp = new double ***[times];     

   for (int t=0;t<times;t++)
   {

         up[t] = new double **[s_level];  
         vp[t] = new double **[s_level];  
         wp[t] = new double **[s_level];      
         for (int k=0;k<s_level;k++)
         {
              up[t][k] = new double *[eta_psi];   
              vp[t][k] = new double *[eta_psi]; 
              wp[t][k] = new double *[eta_psi]; 
              for (int j=0;j<eta_psi;j++)
              {
                   up[t][k][j] = new double [xi_psi];       
                   vp[t][k][j] = new double [xi_psi];   
                   wp[t][k][j] = new double [xi_psi];      
              } 
         }
 
    }



  Variables particle;
  double auxx_s_w=s_w;
  double auxx_s_level=s_level;
  double aux_level=((auxx_s_w-1)/(auxx_s_level-1));


   
   for (int t=0;t<times-1;t++)
   {
       for (int j=0;j<eta_psi-1;j++)
       {
         for (int i=0;i<xi_psi-1;i++)
         {
            double p_s_aux=0;
            for (int k=0;k<s_level-1;k++)
            {
              particle.xi=i;
              particle.eta=j;
              particle.s=p_s_aux;
              particle.t=t;

              UVW(particle);

              up[t][k][j][i]=particle.u;
              vp[t][k][j][i]=particle.v;
              wp[t][k][j][i]=particle.w;

            p_s_aux=p_s_aux+aux_level;
            }
       } 
     }
   }










   for (int t=0;t<times-1;t++)
   {
       for (int j=0;j<eta_psi-1;j++)
       {
         for (int i=0;i<xi_psi-1;i++)
         {
              int k=s_level-1;

              up[t][k][j][i]=0;
              vp[t][k][j][i]=0;
              wp[t][k][j][i]=0;


         }
       } 
   }



  // delete[] u;
  // delete[] v;
  // delete[] w;

}












bool Boundaries(Variables &point)  //xi,eta,s,t, -> Boundaries Domain
{

   bool bounded;

   //psi points
   int xi_i=point.xi;
   int eta_i=point.eta;
   int Q11,Q21,Q12,Q22;
   int bound;


   bounded=false;
   if(point.xi>=0 && point.xi<xi-1 && point.eta>=0 && point.eta<eta-1 && point.t>=0 && point.t<times-1 && point.s>=0 && point.s<s-1)
   {


        Q11=mask_psi[eta_i][xi_i];
        Q21=mask_psi[eta_i][xi_i+1];
        Q12=mask_psi[eta_i+1][xi_i];
        Q22=mask_psi[eta_i+1][xi_i+1];
        bound=Q11*Q12*Q21*Q22;
        if(bound==1){bounded=true;}


       //  bounded=true;
   }
   return bounded;
}














int UVW(Variables &point)  //xi,eta,s,t -> u(xi,eta,s,t), v(xi,eta,s,t), w(xi,eta,s,t)
{




   Hexaedron HexaedronX,HexaedronY,HexaedronZ;

   Hexaedron Hex1X,Hex1Y,Hex1Z;
   Hexaedron Hex2X,Hex2Y,Hex2Z;
   Hexaedron Hex3X,Hex3Y,Hex3Z;



   int xi_i_rho,eta_i_rho,s_i_w;
   double x_rho,y_rho,sg_w;
   double p_xi_rho,p_eta_rho,p_s_w;
   double dt;
   double Q0,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Qt1,Qt2;


  double pdxi,pdeta,pds;
  double J0,J1,J2,J3,J4,J5,J6,J7,J8;
  double M00,M10,M20;
  double M01,M11,M21;
  double M02,M12,M22;

  double J01,J02,J03,J04,J05,J06,J07;
  double U_medio,V_medio,W_medio;
  double det;
  double U0,U1,V0,V1,W0,W1;
  double u0,u1,v0,v1,w0,w1;


   //psi points
   int xi_i=point.xi;
   int eta_i=point.eta;
   int s_i=point.s;
   int t_i=point.t;

   double x=point.xi-xi_i;
   double y=point.eta-eta_i;
   double sg=point.s-s_i;
   double tg=point.t-t_i;

   dxi=x;
   deta=y;
   ds=sg;
   dt=tg;


  //VELOCITIES



  Qt1=u[t_i][s_i][eta_i+1][xi_i];
  Qt2=u[t_i+1][s_i][eta_i+1][xi_i];
  u0=Qt1+((Qt2-Qt1)*(dt));

  Qt1=u[t_i][s_i][eta_i+1][xi_i+1];
  Qt2=u[t_i+1][s_i][eta_i+1][xi_i+1];
  u1=Qt1+((Qt2-Qt1)*(dt));


  Qt1=v[t_i][s_i][eta_i][xi_i+1];
  Qt2=v[t_i+1][s_i][eta_i][xi_i+1];
  v0=Qt1+((Qt2-Qt1)*(dt));

  Qt1=v[t_i][s_i][eta_i+1][xi_i+1];
  Qt2=v[t_i+1][s_i][eta_i+1][xi_i+1];
  v1=Qt1+((Qt2-Qt1)*(dt));


  Qt1=w[t_i][s_i][eta_i+1][xi_i+1];
  Qt2=w[t_i+1][s_i][eta_i+1][xi_i+1];
  w0=Qt1+((Qt2-Qt1)*(dt));

  Qt1=w[t_i][s_i+1][eta_i+1][xi_i+1];
  Qt2=w[t_i+1][s_i+1][eta_i+1][xi_i+1];
  w1=Qt1+((Qt2-Qt1)*(dt));


  //LON,LAT

  Q0=0;
  Q1=relative_lon_psi_Q1_meters[eta_i][xi_i];
  Q3=relative_lon_psi_Q3_meters[eta_i][xi_i];
  Q2=relative_lon_psi_Q2_meters[eta_i][xi_i];

  Q4=Q0;
  Q5=Q1;
  Q6=Q2;
  Q7=Q3;

  HexaedronX.cero=Q0;
  HexaedronX.uno=Q1;
  HexaedronX.dos=Q2;
  HexaedronX.tres=Q3;
  HexaedronX.cuatro=Q4;
  HexaedronX.cinco=Q5;
  HexaedronX.seis=Q6;
  HexaedronX.siete=Q7;

  Q0=0;
  Q1=relative_lat_psi_Q1_meters[eta_i][xi_i];
  Q3=relative_lat_psi_Q3_meters[eta_i][xi_i];
  Q2=relative_lat_psi_Q2_meters[eta_i][xi_i];

  Q4=Q0;
  Q5=Q1;
  Q6=Q2;
  Q7=Q3;

  HexaedronY.cero=Q0;
  HexaedronY.uno=Q1;
  HexaedronY.dos=Q2;
  HexaedronY.tres=Q3;
  HexaedronY.cuatro=Q4;
  HexaedronY.cinco=Q5;
  HexaedronY.seis=Q6;
  HexaedronY.siete=Q7;
  


  //Depth

  Qt1=depth_psi_w[t_i][s_i][eta_i][xi_i];
  Qt2=depth_psi_w[t_i+1][s_i][eta_i][xi_i];
  Q0=Qt1+((Qt2-Qt1)*(dt));

  Qt1=depth_psi_w[t_i][s_i][eta_i][xi_i+1];
  Qt2=depth_psi_w[t_i+1][s_i][eta_i][xi_i+1];
  Q1=Qt1+((Qt2-Qt1)*(dt));

  Qt1=depth_psi_w[t_i][s_i][eta_i+1][xi_i];
  Qt2=depth_psi_w[t_i+1][s_i][eta_i+1][xi_i];
  Q3=Qt1+((Qt2-Qt1)*(dt));

  Qt1=depth_psi_w[t_i][s_i][eta_i+1][xi_i+1];
  Qt2=depth_psi_w[t_i+1][s_i][eta_i+1][xi_i+1];
  Q2=Qt1+((Qt2-Qt1)*(dt));



  Qt1=depth_psi_w[t_i][s_i+1][eta_i][xi_i];
  Qt2=depth_psi_w[t_i+1][s_i+1][eta_i][xi_i];
  Q4=Qt1+((Qt2-Qt1)*(dt));

  Qt1=depth_psi_w[t_i][s_i+1][eta_i][xi_i+1];
  Qt2=depth_psi_w[t_i+1][s_i+1][eta_i][xi_i+1];
  Q5=Qt1+((Qt2-Qt1)*(dt));

  Qt1=depth_psi_w[t_i][s_i+1][eta_i+1][xi_i];
  Qt2=depth_psi_w[t_i+1][s_i+1][eta_i+1][xi_i];
  Q7=Qt1+((Qt2-Qt1)*(dt));

  Qt1=depth_psi_w[t_i][s_i+1][eta_i+1][xi_i+1];
  Qt2=depth_psi_w[t_i+1][s_i+1][eta_i+1][xi_i+1];
  Q6=Qt1+((Qt2-Qt1)*(dt));

  HexaedronZ.cero=Q0;
  HexaedronZ.uno=Q1;
  HexaedronZ.dos=Q2;
  HexaedronZ.tres=Q3;
  HexaedronZ.cuatro=Q4;
  HexaedronZ.cinco=Q5;
  HexaedronZ.seis=Q6;
  HexaedronZ.siete=Q7;




  //HEXAEDRONES 1/2


  //HEXADRON_1

  Hex1X.cero=HexaedronX.cero;
  Hex1X.tres=HexaedronX.tres;
  Hex1X.cuatro=HexaedronX.cuatro;
  Hex1X.siete=HexaedronX.siete;

  Hex1X.uno=(HexaedronX.cero+HexaedronX.uno)*0.5;
  Hex1X.dos=(HexaedronX.tres+HexaedronX.dos)*0.5;
  Hex1X.seis=(HexaedronX.siete+HexaedronX.seis)*0.5;
  Hex1X.cinco=(HexaedronX.cuatro+HexaedronX.cinco)*0.5;

  Hex1Y.cero=HexaedronY.cero;
  Hex1Y.tres=HexaedronY.tres;
  Hex1Y.cuatro=HexaedronY.cuatro;
  Hex1Y.siete=HexaedronY.siete;

  Hex1Y.uno=(HexaedronY.cero+HexaedronY.uno)*0.5;
  Hex1Y.dos=(HexaedronY.tres+HexaedronY.dos)*0.5;
  Hex1Y.seis=(HexaedronY.siete+HexaedronY.seis)*0.5;
  Hex1Y.cinco=(HexaedronY.cuatro+HexaedronY.cinco)*0.5;

  Hex1Z.cero=HexaedronZ.cero;
  Hex1Z.tres=HexaedronZ.tres;
  Hex1Z.cuatro=HexaedronZ.cuatro;
  Hex1Z.siete=HexaedronZ.siete;

  Hex1Z.uno=(HexaedronZ.cero+HexaedronZ.uno)*0.5;
  Hex1Z.dos=(HexaedronZ.tres+HexaedronZ.dos)*0.5;
  Hex1Z.seis=(HexaedronZ.siete+HexaedronZ.seis)*0.5;
  Hex1Z.cinco=(HexaedronZ.cuatro+HexaedronZ.cinco)*0.5;


  //HEXADRON_2

  Hex2X.cero=HexaedronX.cero;
  Hex2X.uno=HexaedronX.uno;
  Hex2X.cuatro=HexaedronX.cuatro;
  Hex2X.cinco=HexaedronX.cinco;

  Hex2X.tres=(HexaedronX.cero+HexaedronX.tres)*0.5;
  Hex2X.dos=(HexaedronX.uno+HexaedronX.dos)*0.5;
  Hex2X.seis=(HexaedronX.cinco+HexaedronX.seis)*0.5;
  Hex2X.siete=(HexaedronX.cuatro+HexaedronX.siete)*0.5;


  Hex2Y.cero=HexaedronY.cero;
  Hex2Y.uno=HexaedronY.uno;
  Hex2Y.cuatro=HexaedronY.cuatro;
  Hex2Y.cinco=HexaedronY.cinco;

  Hex2Y.tres=(HexaedronY.cero+HexaedronY.tres)*0.5;
  Hex2Y.dos=(HexaedronY.uno+HexaedronY.dos)*0.5;
  Hex2Y.seis=(HexaedronY.cinco+HexaedronY.seis)*0.5;
  Hex2Y.siete=(HexaedronY.cuatro+HexaedronY.siete)*0.5;


  Hex2Z.cero=HexaedronZ.cero;
  Hex2Z.uno=HexaedronZ.uno;
  Hex2Z.cuatro=HexaedronZ.cuatro;
  Hex2Z.cinco=HexaedronZ.cinco;

  Hex2Z.tres=(HexaedronZ.cero+HexaedronZ.tres)*0.5;
  Hex2Z.dos=(HexaedronZ.uno+HexaedronZ.dos)*0.5;
  Hex2Z.seis=(HexaedronZ.cinco+HexaedronZ.seis)*0.5;
  Hex2Z.siete=(HexaedronZ.cuatro+HexaedronZ.siete)*0.5;



  //HEXADRON_3


  Hex3X.cero=HexaedronX.cero;
  Hex3X.uno=HexaedronX.uno;
  Hex3X.dos=HexaedronX.dos;
  Hex3X.tres=HexaedronX.tres;

  Hex3X.cuatro=(HexaedronX.cero+HexaedronX.cuatro)*0.5;
  Hex3X.cinco=(HexaedronX.uno+HexaedronX.cinco)*0.5;
  Hex3X.seis=(HexaedronX.dos+HexaedronX.seis)*0.5;
  Hex3X.siete=(HexaedronX.tres+HexaedronX.siete)*0.5;


  Hex3Y.cero=HexaedronY.cero;
  Hex3Y.uno=HexaedronY.uno;
  Hex3Y.dos=HexaedronY.dos;
  Hex3Y.tres=HexaedronY.tres;

  Hex3Y.cuatro=(HexaedronY.cero+HexaedronY.cuatro)*0.5;
  Hex3Y.cinco=(HexaedronY.uno+HexaedronY.cinco)*0.5;
  Hex3Y.seis=(HexaedronY.dos+HexaedronY.seis)*0.5;
  Hex3Y.siete=(HexaedronY.tres+HexaedronY.siete)*0.5;


  Hex3Z.cero=HexaedronZ.cero;
  Hex3Z.uno=HexaedronZ.uno;
  Hex3Z.dos=HexaedronZ.dos;
  Hex3Z.tres=HexaedronZ.tres;

  Hex3Z.cuatro=(HexaedronZ.cero+HexaedronZ.cuatro)*0.5;
  Hex3Z.cinco=(HexaedronZ.uno+HexaedronZ.cinco)*0.5;
  Hex3Z.seis=(HexaedronZ.dos+HexaedronZ.seis)*0.5;
  Hex3Z.siete=(HexaedronZ.tres+HexaedronZ.siete)*0.5;


  //FIN HEXAEDRONES 1/2



  //JACOBIAN 3D
  pdxi=dxi;
  pdeta=deta;
  pds=ds;

  Jacobian(HexaedronX,HexaedronY,HexaedronZ,pdxi,pdeta,pds);

  det=(J[0][0]*J[1][1]*J[2][2])+(J[1][0]*J[2][1]*J[0][2])+(J[0][1]*J[1][2]*J[2][0])-(J[2][0]*J[1][1]*J[0][2])-(J[1][0]*J[0][1]*J[2][2])-(J[2][1]*J[1][2]*J[0][0]);


  //J0
  pdxi=0;
  pdeta=deta;
  pds=ds;

  Jacobian(HexaedronX,HexaedronY,HexaedronZ,pdxi,pdeta,pds);

  M00=J[1][1]*J[2][2]-J[2][1]*J[1][2];
  M10=J[0][1]*J[2][2]-J[2][1]*J[0][2];
  M20=J[0][1]*J[1][2]-J[1][1]*J[0][2];

  J0=sqrt(M00*M00+M10*M10+M20*M20);

  //J1
  pdxi=1;
  pdeta=deta;
  pds=ds;

  Jacobian(HexaedronX,HexaedronY,HexaedronZ,pdxi,pdeta,pds);

  M00=J[1][1]*J[2][2]-J[2][1]*J[1][2];
  M10=J[0][1]*J[2][2]-J[2][1]*J[0][2];
  M20=J[0][1]*J[1][2]-J[1][1]*J[0][2];

  J1=sqrt(M00*M00+M10*M10+M20*M20);

  //J2
  pdxi=0.5;
  pdeta=deta;
  pds=ds;

  Jacobian(HexaedronX,HexaedronY,HexaedronZ,pdxi,pdeta,pds);

  M00=J[1][1]*J[2][2]-J[2][1]*J[1][2];
  M10=J[0][1]*J[2][2]-J[2][1]*J[0][2];
  M20=J[0][1]*J[1][2]-J[1][1]*J[0][2];

  J2=sqrt(M00*M00+M10*M10+M20*M20);


  //J3
  pdxi=dxi;
  pdeta=0;
  pds=ds;

  Jacobian(HexaedronX,HexaedronY,HexaedronZ,pdxi,pdeta,pds);

  M01=J[1][0]*J[2][2]-J[2][0]*J[1][2];
  M11=J[0][0]*J[2][2]-J[2][0]*J[0][2];
  M21=J[0][0]*J[1][2]-J[1][0]*J[0][2];

  J3=sqrt(M01*M01+M11*M11+M21*M21);

  //J4
  pdxi=dxi;
  pdeta=1;
  pds=ds;

  Jacobian(HexaedronX,HexaedronY,HexaedronZ,pdxi,pdeta,pds);

  M01=J[1][0]*J[2][2]-J[2][0]*J[1][2];
  M11=J[0][0]*J[2][2]-J[2][0]*J[0][2];
  M21=J[0][0]*J[1][2]-J[1][0]*J[0][2];

  J4=sqrt(M01*M01+M11*M11+M21*M21);

  //J5
  pdxi=dxi;
  pdeta=0.5;
  pds=ds;

  Jacobian(HexaedronX,HexaedronY,HexaedronZ,pdxi,pdeta,pds);

  M01=J[1][0]*J[2][2]-J[2][0]*J[1][2];
  M11=J[0][0]*J[2][2]-J[2][0]*J[0][2];
  M21=J[0][0]*J[1][2]-J[1][0]*J[0][2];

  J5=sqrt(M01*M01+M11*M11+M21*M21);


  //J6
  pdxi=dxi;
  pdeta=deta;
  pds=0;

  Jacobian(HexaedronX,HexaedronY,HexaedronZ,pdxi,pdeta,pds);

  M02=J[1][0]*J[2][1]-J[2][0]*J[1][1];
  M12=J[0][0]*J[2][1]-J[2][0]*J[0][1];
  M22=J[0][0]*J[1][1]-J[1][0]*J[0][1];

  J6=sqrt(M02*M02+M12*M12+M22*M22);


  //J7
  pdxi=dxi;
  pdeta=deta;
  pds=1;

  Jacobian(HexaedronX,HexaedronY,HexaedronZ,pdxi,pdeta,pds);

  M02=J[1][0]*J[2][1]-J[2][0]*J[1][1];
  M12=J[0][0]*J[2][1]-J[2][0]*J[0][1];
  M22=J[0][0]*J[1][1]-J[1][0]*J[0][1];

  J7=sqrt(M02*M02+M12*M12+M22*M22);


  //J8
  pdxi=dxi;
  pdeta=deta;
  pds=0.5;

  Jacobian(HexaedronX,HexaedronY,HexaedronZ,pdxi,pdeta,pds);

  M02=J[1][0]*J[2][1]-J[2][0]*J[1][1];
  M12=J[0][0]*J[2][1]-J[2][0]*J[0][1];
  M22=J[0][0]*J[1][1]-J[1][0]*J[0][1];

  J8=sqrt(M02*M02+M12*M12+M22*M22);





  //Jacobianos 1/2  -> U_1/2

  


  //J_1/2 J01


  pdxi=0;
  pdeta=0.5;
  pds=0.5;

  Jacobian(Hex1X,Hex1Y,Hex1Z,pdxi,pdeta,pds);


  M00=J[1][1]*J[2][2]-J[2][1]*J[1][2];
  M10=J[0][1]*J[2][2]-J[2][1]*J[0][2];
  M20=J[0][1]*J[1][2]-J[1][1]*J[0][2];

  J01=sqrt(M00*M00+M10*M10+M20*M20);


  //J_1/2  J02


  pdxi=0.5;
  pdeta=0;
  pds=0.5;

  Jacobian(Hex1X,Hex1Y,Hex1Z,pdxi,pdeta,pds);


  M01=J[1][0]*J[2][2]-J[2][0]*J[1][2];
  M11=J[0][0]*J[2][2]-J[2][0]*J[0][2];
  M21=J[0][0]*J[1][2]-J[1][0]*J[0][2];

  J02=sqrt(M01*M01+M11*M11+M21*M21);

  //J_1/2  J03


  pdxi=0.5;
  pdeta=1;
  pds=0.5;

  Jacobian(Hex1X,Hex1Y,Hex1Z,pdxi,pdeta,pds);

  M01=J[1][0]*J[2][2]-J[2][0]*J[1][2];
  M11=J[0][0]*J[2][2]-J[2][0]*J[0][2];
  M21=J[0][0]*J[1][2]-J[1][0]*J[0][2];

  J03=sqrt(M01*M01+M11*M11+M21*M21);

  //J_1/2  J04


  pdxi=0.5;
  pdeta=0.5;
  pds=0;

  Jacobian(Hex1X,Hex1Y,Hex1Z,pdxi,pdeta,pds);

  M02=J[1][0]*J[2][1]-J[2][0]*J[1][1];
  M12=J[0][0]*J[2][1]-J[2][0]*J[0][1];
  M22=J[0][0]*J[1][1]-J[1][0]*J[0][1];

  J04=sqrt(M02*M02+M12*M12+M22*M22);


  //J_1/2  J05


  pdxi=0.5;
  pdeta=0.5;
  pds=1;

  Jacobian(Hex1X,Hex1Y,Hex1Z,pdxi,pdeta,pds);

  M02=J[1][0]*J[2][1]-J[2][0]*J[1][1];
  M12=J[0][0]*J[2][1]-J[2][0]*J[0][1];
  M22=J[0][0]*J[1][1]-J[1][0]*J[0][1];

  J05=sqrt(M02*M02+M12*M12+M22*M22);

  //J_1/2  J06


  pdxi=1;
  pdeta=deta;
  pds=ds;

  Jacobian(Hex1X,Hex1Y,Hex1Z,pdxi,pdeta,pds);

  M00=J[1][1]*J[2][2]-J[2][1]*J[1][2];
  M10=J[0][1]*J[2][2]-J[2][1]*J[0][2];
  M20=J[0][1]*J[1][2]-J[1][1]*J[0][2];

  J06=sqrt(M00*M00+M10*M10+M20*M20);


  //J_1/2  J07


  pdxi=1;
  pdeta=0.5;
  pds=0.5;

  Jacobian(Hex1X,Hex1Y,Hex1Z,pdxi,pdeta,pds);

  M00=J[1][1]*J[2][2]-J[2][1]*J[1][2];
  M10=J[0][1]*J[2][2]-J[2][1]*J[0][2];
  M20=J[0][1]*J[1][2]-J[1][1]*J[0][2];

  J07=sqrt(M00*M00+M10*M10+M20*M20);



  U_medio=(u0*J01+v0*J02-v1*J03+w0*J04-w1*J05);
  U_medio=U_medio*(J06/J07);







  //Jacobianos 1/2  -> V_1/2

  


  //J_1/2 J01


  pdxi=0;
  pdeta=0.5;
  pds=0.5;

  Jacobian(Hex2X,Hex2Y,Hex2Z,pdxi,pdeta,pds);

  M00=J[1][1]*J[2][2]-J[2][1]*J[1][2];
  M10=J[0][1]*J[2][2]-J[2][1]*J[0][2];
  M20=J[0][1]*J[1][2]-J[1][1]*J[0][2];

  J01=sqrt(M00*M00+M10*M10+M20*M20);



  //J_1/2 J02


  pdxi=1;
  pdeta=0.5;
  pds=0.5;

  Jacobian(Hex2X,Hex2Y,Hex2Z,pdxi,pdeta,pds);

  M00=J[1][1]*J[2][2]-J[2][1]*J[1][2];
  M10=J[0][1]*J[2][2]-J[2][1]*J[0][2];
  M20=J[0][1]*J[1][2]-J[1][1]*J[0][2];

  J02=sqrt(M00*M00+M10*M10+M20*M20);


  //J_1/2 J03


  pdxi=0.5;
  pdeta=0;
  pds=0.5;

  Jacobian(Hex2X,Hex2Y,Hex2Z,pdxi,pdeta,pds);

  M01=J[1][0]*J[2][2]-J[2][0]*J[1][2];
  M11=J[0][0]*J[2][2]-J[2][0]*J[0][2];
  M21=J[0][0]*J[1][2]-J[1][0]*J[0][2];

  J03=sqrt(M01*M01+M11*M11+M21*M21);

  //J_1/2 J04


  pdxi=0.5;
  pdeta=0.5;
  pds=0;

  Jacobian(Hex2X,Hex2Y,Hex2Z,pdxi,pdeta,pds);

  M02=J[1][0]*J[2][1]-J[2][0]*J[1][1];
  M12=J[0][0]*J[2][1]-J[2][0]*J[0][1];
  M22=J[0][0]*J[1][1]-J[1][0]*J[0][1];

  J04=sqrt(M02*M02+M12*M12+M22*M22);


  //J_1/2 J05


  pdxi=0.5;
  pdeta=0.5;
  pds=1;

  Jacobian(Hex2X,Hex2Y,Hex2Z,pdxi,pdeta,pds);

  M02=J[1][0]*J[2][1]-J[2][0]*J[1][1];
  M12=J[0][0]*J[2][1]-J[2][0]*J[0][1];
  M22=J[0][0]*J[1][1]-J[1][0]*J[0][1];

  J05=sqrt(M02*M02+M12*M12+M22*M22);


  //J_1/2 J06


  pdxi=dxi;
  pdeta=1;
  pds=ds;

  Jacobian(Hex2X,Hex2Y,Hex2Z,pdxi,pdeta,pds);

  M01=J[1][0]*J[2][2]-J[2][0]*J[1][2];
  M11=J[0][0]*J[2][2]-J[2][0]*J[0][2];
  M21=J[0][0]*J[1][2]-J[1][0]*J[0][2];

  J06=sqrt(M01*M01+M11*M11+M21*M21);

  //J_1/2 J07


  pdxi=0.5;
  pdeta=1;
  pds=0.5;

  Jacobian(Hex2X,Hex2Y,Hex2Z,pdxi,pdeta,pds);

  M01=J[1][0]*J[2][2]-J[2][0]*J[1][2];
  M11=J[0][0]*J[2][2]-J[2][0]*J[0][2];
  M21=J[0][0]*J[1][2]-J[1][0]*J[0][2];

  J07=sqrt(M01*M01+M11*M11+M21*M21);



  V_medio=(u0*J01-u1*J02+v0*J03+w0*J04-w1*J05);
  V_medio=V_medio*(J06/J07);


  //Jacobianos 1/2  -> W_1/2

  


  //J_1/2 J01


  pdxi=0;
  pdeta=0.5;
  pds=0.5;

  Jacobian(Hex3X,Hex3Y,Hex3Z,pdxi,pdeta,pds);

  M00=J[1][1]*J[2][2]-J[2][1]*J[1][2];
  M10=J[0][1]*J[2][2]-J[2][1]*J[0][2];
  M20=J[0][1]*J[1][2]-J[1][1]*J[0][2];

  J01=sqrt(M00*M00+M10*M10+M20*M20);


  //J_1/2 J02


  pdxi=1;
  pdeta=0.5;
  pds=0.5;

  Jacobian(Hex3X,Hex3Y,Hex3Z,pdxi,pdeta,pds);

  M00=J[1][1]*J[2][2]-J[2][1]*J[1][2];
  M10=J[0][1]*J[2][2]-J[2][1]*J[0][2];
  M20=J[0][1]*J[1][2]-J[1][1]*J[0][2];

  J02=sqrt(M00*M00+M10*M10+M20*M20);

  //J_1/2 J03


  pdxi=0.5;
  pdeta=0;
  pds=0.5;

  Jacobian(Hex3X,Hex3Y,Hex3Z,pdxi,pdeta,pds);

  M01=J[1][0]*J[2][2]-J[2][0]*J[1][2];
  M11=J[0][0]*J[2][2]-J[2][0]*J[0][2];
  M21=J[0][0]*J[1][2]-J[1][0]*J[0][2];

  J03=sqrt(M01*M01+M11*M11+M21*M21);

  //J_1/2 J04


  pdxi=0.5;
  pdeta=1;
  pds=0.5;

  Jacobian(Hex3X,Hex3Y,Hex3Z,pdxi,pdeta,pds);

  M01=J[1][0]*J[2][2]-J[2][0]*J[1][2];
  M11=J[0][0]*J[2][2]-J[2][0]*J[0][2];
  M21=J[0][0]*J[1][2]-J[1][0]*J[0][2];

  J04=sqrt(M01*M01+M11*M11+M21*M21);

  //J_1/2 J05


  pdxi=0.5;
  pdeta=0.5;
  pds=0;

  Jacobian(Hex3X,Hex3Y,Hex3Z,pdxi,pdeta,pds);

  M02=J[1][0]*J[2][1]-J[2][0]*J[1][1];
  M12=J[0][0]*J[2][1]-J[2][0]*J[0][1];
  M22=J[0][0]*J[1][1]-J[1][0]*J[0][1];

  J05=sqrt(M02*M02+M12*M12+M22*M22);


  //J_1/2 J06


  pdxi=dxi;
  pdeta=deta;
  pds=1;

  Jacobian(Hex3X,Hex3Y,Hex3Z,pdxi,pdeta,pds);

  M02=J[1][0]*J[2][1]-J[2][0]*J[1][1];
  M12=J[0][0]*J[2][1]-J[2][0]*J[0][1];
  M22=J[0][0]*J[1][1]-J[1][0]*J[0][1];

  J06=sqrt(M02*M02+M12*M12+M22*M22);


  //J_1/2 J07


  pdxi=0.5;
  pdeta=0.5;
  pds=1;

  Jacobian(Hex3X,Hex3Y,Hex3Z,pdxi,pdeta,pds);

  M02=J[1][0]*J[2][1]-J[2][0]*J[1][1];
  M12=J[0][0]*J[2][1]-J[2][0]*J[0][1];
  M22=J[0][0]*J[1][1]-J[1][0]*J[0][1];

  J07=sqrt(M02*M02+M12*M12+M22*M22);



  W_medio=(u0*J01-u1*J02+v0*J03-v1*J04+w0*J05);
  W_medio=W_medio*(J06/J07);


 //FLUXES

  U0=u0*J0;
  U1=u1*J1;
  V0=v0*J3;
  V1=v1*J4;
  W0=w0*J6;
  W1=w1*J7;
  

 //VELOCITIES

  point.u=((((2*dxi*dxi)-(3*dxi)+1)*U0)+(((-4*dxi*dxi)+(4*dxi))*U_medio)+(((2*dxi*dxi)-(dxi))*U1))/det;
  point.v=((((2*deta*deta)-(3*deta)+1)*V0)+(((-4*deta*deta)+(4*deta))*V_medio)+(((2*deta*deta)-(deta))*V1))/det;
  point.w=((((2*ds*ds)-(3*ds)+1)*W0)+(((-4*ds*ds)+(4*ds))*W_medio)+(((2*ds*ds)-(ds))*W1))/det;




  point.u=point.u*secs_per_unit_time;
  point.v=point.v*secs_per_unit_time;
  point.w=point.w*secs_per_unit_time;


 

}








int Jacobian(Hexaedron &HX,Hexaedron &HY,Hexaedron &HZ,double &dxi,double &deta,double &ds)
{

  double dx_dxi,dx_deta,dx_ds;
  double dy_dxi,dy_deta,dy_ds;
  double dz_dxi,dz_deta,dz_ds;


  dx_dxi=(-(1-deta)*(1-ds)*HX.cero)+((1-deta)*(1-ds)*HX.uno)+(deta*(1-ds)*HX.dos)+(-deta*(1-ds)*HX.tres)+(-(1-deta)*(ds)*HX.cuatro)+((1-deta)*(ds)*HX.cinco)+((deta)*(ds)*HX.seis)+(-(deta)*(ds)*HX.siete);

  dx_deta=(-(1-dxi)*(1-ds)*HX.cero)+(-(dxi)*(1-ds)*HX.uno)+((dxi)*(1-ds)*HX.dos)+((1-dxi)*(1-ds)*HX.tres)+(-(1-dxi)*(ds)*HX.cuatro)+(-(dxi)*(ds)*HX.cinco)+((dxi)*(ds)*HX.seis)+((1-dxi)*(ds)*HX.siete);

  dx_ds=(-(1-dxi)*(1-deta)*HX.cero)+(-(dxi)*(1-deta)*HX.uno)+(-(dxi)*(deta)*HX.dos)+(-(1-dxi)*(deta)*HX.tres)+((1-dxi)*(1-deta)*HX.cuatro)+((dxi)*(1-deta)*HX.cinco)+((dxi)*(deta)*HX.seis)+((1-dxi)*(deta)*HX.siete);



  dy_dxi=(-(1-deta)*(1-ds)*HY.cero)+((1-deta)*(1-ds)*HY.uno)+(deta*(1-ds)*HY.dos)+(-deta*(1-ds)*HY.tres)+(-(1-deta)*(ds)*HY.cuatro)+((1-deta)*(ds)*HY.cinco)+((deta)*(ds)*HY.seis)+(-(deta)*(ds)*HY.siete);

  dy_deta=(-(1-dxi)*(1-ds)*HY.cero)+(-(dxi)*(1-ds)*HY.uno)+((dxi)*(1-ds)*HY.dos)+((1-dxi)*(1-ds)*HY.tres)+(-(1-dxi)*(ds)*HY.cuatro)+(-(dxi)*(ds)*HY.cinco)+((dxi)*(ds)*HY.seis)+((1-dxi)*(ds)*HY.siete);

  dy_ds=(-(1-dxi)*(1-deta)*HY.cero)+(-(dxi)*(1-deta)*HY.uno)+(-(dxi)*(deta)*HY.dos)+(-(1-dxi)*(deta)*HY.tres)+((1-dxi)*(1-deta)*HY.cuatro)+((dxi)*(1-deta)*HY.cinco)+((dxi)*(deta)*HY.seis)+((1-dxi)*(deta)*HY.siete);



  dz_dxi=(-(1-deta)*(1-ds)*HZ.cero)+((1-deta)*(1-ds)*HZ.uno)+(deta*(1-ds)*HZ.dos)+(-deta*(1-ds)*HZ.tres)+(-(1-deta)*(ds)*HZ.cuatro)+((1-deta)*(ds)*HZ.cinco)+((deta)*(ds)*HZ.seis)+(-(deta)*(ds)*HZ.siete);

  dz_deta=(-(1-dxi)*(1-ds)*HZ.cero)+(-(dxi)*(1-ds)*HZ.uno)+((dxi)*(1-ds)*HZ.dos)+((1-dxi)*(1-ds)*HZ.tres)+(-(1-dxi)*(ds)*HZ.cuatro)+(-(dxi)*(ds)*HZ.cinco)+((dxi)*(ds)*HZ.seis)+((1-dxi)*(ds)*HZ.siete);

  dz_ds=(-(1-dxi)*(1-deta)*HZ.cero)+(-(dxi)*(1-deta)*HZ.uno)+(-(dxi)*(deta)*HZ.dos)+(-(1-dxi)*(deta)*HZ.tres)+((1-dxi)*(1-deta)*HZ.cuatro)+((dxi)*(1-deta)*HZ.cinco)+((dxi)*(deta)*HZ.seis)+((1-dxi)*(deta)*HZ.siete);


  J[0][0]=dx_dxi;
  J[0][1]=dx_deta;
  J[0][2]=dx_ds;

  J[1][0]=dy_dxi;
  J[1][1]=dy_deta;
  J[1][2]=dy_ds;

  J[2][0]=dz_dxi;
  J[2][1]=dz_deta;
  J[2][2]=dz_ds;



}












int conversion_uv()
{






  mask_psi = new double *[eta_psi];
   for (int i = 0; i < eta_psi; i++)
   {
     mask_psi[i] = new double [xi_psi];
  
   }




          double Q11,Q12,Q21,Q22;

          for (int j = 0; j < eta_psi; j++)
          {
             for (int i = 0; i < xi_psi; i++)
             {

                Q11=mask_rho[j][i];
                Q21=mask_rho[j][i+1];
                Q12=mask_rho[j+1][i];
                Q22=mask_rho[j+1][i+1];
                mask_psi[j][i]=Q11*Q12*Q21*Q22;

             }
          }





   for(int tk=0;tk<times;tk++)
   {
     for(int k=0;k<s_w;k++)
     {
       for(int j=0;j<eta_rho;j++)
       {
         for(int i=0;i<xi_rho;i++)
         {
             w[tk][k][j][i]=(w[tk][k][j][i]+sinking_velocity); //*mask_rho[j][i];
         }
       } 
     }
   }









  relative_lon_psi_Q1_meters = new double *[eta_psi];
  relative_lon_psi_Q3_meters = new double *[eta_psi];
  relative_lon_psi_Q2_meters = new double *[eta_psi];
  relative_lat_psi_Q1_meters = new double *[eta_psi];
  relative_lat_psi_Q3_meters = new double *[eta_psi];
  relative_lat_psi_Q2_meters = new double *[eta_psi];

   for (int i = 0; i < eta_psi; i++)
   {
     relative_lon_psi_Q1_meters[i] = new double [xi_psi];
     relative_lon_psi_Q3_meters[i] = new double [xi_psi];
     relative_lon_psi_Q2_meters[i] = new double [xi_psi];
     relative_lat_psi_Q1_meters[i] = new double [xi_psi];
     relative_lat_psi_Q3_meters[i] = new double [xi_psi];
     relative_lat_psi_Q2_meters[i] = new double [xi_psi];
  
   }




     double dlongitude,dlatitude,longitude3,latitude3,longitude2,latitude2,longitude1,latitude1,longitude0,latitude0;

 
       for(int j=0;j<eta_psi-1;j++)
       {
         for(int i=0;i<xi_psi-1;i++)
         {

             longitude0=lon_psi[j][i];
             latitude0=lat_psi[j][i];

             longitude1=lon_psi[j][i+1];
             latitude1=lat_psi[j][i+1];

             dlongitude=distance(longitude0,latitude1,longitude1,latitude1);
             dlatitude=distance(longitude1,latitude0,longitude1,latitude1);
             if(latitude1<latitude0){dlatitude=-dlatitude;}
             relative_lon_psi_Q1_meters[j][i]=dlongitude;
             relative_lat_psi_Q1_meters[j][i]=dlatitude;

             longitude2=lon_psi[j+1][i+1];
             latitude2=lat_psi[j+1][i+1];

             dlongitude=distance(longitude0,latitude2,longitude2,latitude2);
             dlatitude=distance(longitude2,latitude0,longitude2,latitude2);
             if(latitude2<latitude0){dlatitude=-dlatitude;}
             relative_lon_psi_Q2_meters[j][i]=dlongitude;
             relative_lat_psi_Q2_meters[j][i]=dlatitude;

             longitude3=lon_psi[j+1][i];
             latitude3=lat_psi[j+1][i];

             dlongitude=distance(longitude0,latitude3,longitude3,latitude3);
             dlatitude=distance(longitude3,latitude0,longitude3,latitude3);
             if(latitude3<latitude0){dlatitude=-dlatitude;}
             relative_lon_psi_Q3_meters[j][i]=dlongitude;
             relative_lat_psi_Q3_meters[j][i]=dlatitude;


         }
       } 





}




double distance(double &lon1,double &lat1,double &lon2,double &lat2) 
{ 
     double ans;
     

     ans=rearth*2*asin(sqrt((sin(rads*(lat1-lat2)/2))*(sin(rads*(lat1-lat2)/2))+cos(rads*lat1)*cos(rads*lat2)*(sin(rads*(lon2-lon1)/2))*(sin(rads*(lon2-lon1)/2))));
  
    return ans; 
} 








int conversion(double &lon1,double &lat1,double &meters_u,double &meters_v,double &degrees_u,double &degrees_v)
{

  double d2;


  d2=2*asin((sin((meters_u/(2*rearth))))/(cos(rads*lat1)));
  degrees_u=d2*degrees;


  d2=meters_v/rearth;
  degrees_v=d2*degrees;


}






int Scoordinate()
{

  double Csur,Cs_w[s_w],Cs_rho[s_rho],Cbot,mu,Sp;
  double omega_w[s_w],omega_rho[s_rho];
  double nauxw=s_w;
  double naux=s_rho;
  


  depth_rho = new double ***[times];
  depth_w = new double ***[times];
  depth_psi_w = new double ***[times];
   for (int t = 0; t < times; t++)
   {
      depth_rho[t] = new double **[s_rho];
      depth_w[t] = new double **[s_w];
      depth_psi_w[t] = new double **[s_w];
      for (int k = 0; k < s_rho; k++)
      {
         depth_rho[t][k] = new double *[eta_rho];
         for (int j = 0; j < eta_rho; j++)
         {
            depth_rho[t][k][j] = new double [xi_rho];
         }
      }

      for (int k = 0; k < s_w; k++)
      {
         depth_w[t][k] = new double *[eta_rho];   
         depth_psi_w[t][k] = new double *[eta_psi];
         for (int j = 0; j < eta_psi; j++)
         {
            depth_psi_w[t][k][j] = new double [xi_psi];
         }

         for (int j = 0; j < eta_rho; j++)
         {
            depth_w[t][k][j] = new double [xi_rho];
         }
      }
   }

  for(int i=0; i<s_w;i++){ omega_w[i]=(i-naux)/naux;}
  for(int i=0; i<s_rho;i++){ omega_rho[i]=((i+1)-naux-0.5)/(naux);}




  double theta_s,theta_b,hc;

  //Vstretching = 4 (A. Shchepetkin (2010) UCLA-ROMS 
  theta_s=6;
  theta_b=2;
  hc=120;
  mu=0;


  for(int i=0; i<s_w;i++)
  {
   Csur=(1-cosh(omega_w[i]*theta_s))/(cosh(theta_s)-1);
   Cs_w[i]=(exp(theta_b*Csur)-1)/(1-exp(-theta_b));
  }


  for(int i=0; i<s_rho;i++)
  {
   Csur=(1-cosh(omega_rho[i]*theta_s))/(cosh(theta_s)-1);
   Cs_rho[i]=(exp(theta_b*Csur)-1)/(1-exp(-theta_b));
  }


   for (int t = 0; t < times; t++)
   {
      for (int k = 0; k < s_rho; k++)
      {
          for (int j = 0; j < eta_rho; j++)
          {
             for (int i = 0; i < xi_rho; i++)
             {
                 depth_rho[t][k][j][i]=free_surface[t][j][i]+((free_surface[t][j][i]+h[j][i])*(((hc*omega_rho[k])+(h[j][i]*Cs_rho[k]))/(hc+h[j][i])));
             }
          }
       }


      for (int k = 0; k < s_w; k++)
      {
          for (int j = 0; j < eta_rho; j++)
          {
             for (int i = 0; i < xi_rho; i++)
             {
                 depth_w[t][k][j][i]=free_surface[t][j][i]+((free_surface[t][j][i]+h[j][i])*(((hc*omega_w[k])+(h[j][i]*Cs_w[k]))/(hc+h[j][i])));
             }
          }
       }

   }


   double daux;
   double Q11,Q12,Q21,Q22;

   for (int t = 0; t < times; t++)
   {
      for (int k = 0; k < s_w; k++)
      {
          for (int j = 0; j < eta_psi; j++)
          {
             for (int i = 0; i < xi_psi; i++)
             {

                Q11=depth_w[t][k][j][i];
                Q21=depth_w[t][k][j][i+1];
                Q12=depth_w[t][k][j+1][i];
                Q22=depth_w[t][k][j+1][i+1];
                daux=0.5;
                depth_psi_w[t][k][j][i]=(Q11*(1-daux)*(1-daux))+(Q21*(daux)*(1-daux))+(Q12*(1-daux)*(daux))+(Q22*(daux)*(daux));

             }
          }
       }
   }





}












int ReadingVelocityFlow()
{

//Parametros archivos
 initialt=InitialPeriod; 




 int numarchivos;

 double auxiliar=timesim;
 auxiliar=auxiliar/30;
 numarchivos=auxiliar;
 numarchivos=numarchivos+1;
 times=numarchivos*ntime;
 //cout<<"You have a number of: "<<numarchivos*30<<" days available for the simulation"<<endl;  


  NcVar *varb;
  double valor;

  free_surface = new double **[times];   
  u = new double ***[times];  
  v = new double ***[times]; 
  w = new double ***[times];       


   for (int t=0;t<times;t++)
   {
         free_surface[t] = new double *[eta_rho];   
         for (int j=0;j<eta_rho;j++)
         {
              free_surface[t][j] = new double [xi_rho];   
         }

         u[t] = new double **[s_rho];  
         v[t] = new double **[s_rho];       
         for (int k=0;k<s_rho;k++)
         {
              u[t][k] = new double *[eta_rho];   
              v[t][k] = new double *[eta_v];    
              for (int j=0;j<eta_rho;j++)
              {
                   u[t][k][j] = new double [xi_u];       
              } 
              for (int j=0;j<eta_v;j++)
              { 
                   v[t][k][j] = new double [xi_rho];     
              } 
         }
 


         w[t] = new double **[s_w];  
         for (int k=0;k<s_w;k++)
         {
              w[t][k] = new double *[eta_rho];      
              for (int j=0;j<eta_rho;j++)
              {
                   w[t][k][j] = new double [xi_rho];      
              } 
         }
   }






  int i1,i2,i3,i4,i5;
  float *uG, *vG, *wG, *free_surfaceG;
  free_surfaceG = new float [ntime*xi_rho*eta_rho];       
  uG = new float [ntime*xi_u*eta_rho*s_rho];  
  vG = new float [ntime*xi_rho*eta_v*s_rho];  
  wG = new float [ntime*xi_rho*eta_rho*s_w];   





 int t=0;
 int auxi=initialt;
 for(int kaux1=0;kaux1<numarchivos;kaux1++)
 {

     sprintf(ncfile, "/data/geo/escola/roms_canaryislands/cb_2009_3km_42/roms_avg.0%1u.nc",auxi);
     NcFile dataFile(ncfile,NcFile::ReadOnly);
     if(!dataFile.is_valid())
     return NC_ERR;






                      //Read free_surface
                      if (!(varb= dataFile.get_var("zeta")))
	                    return NC_ERR; 
                      if (!varb->set_cur(0,0,0))
	                    return NC_ERR; 
                      if (!varb->get(&free_surfaceG[0],ntime,eta_rho,xi_rho))
	                    return NC_ERR; 
 
                      //Read u
                      if (!(varb= dataFile.get_var("u")))
	                    return NC_ERR; 
                      if (!varb->set_cur(0,0,0,0))
	                    return NC_ERR; 
                      if (!varb->get(&uG[0],ntime,s_rho,eta_rho,xi_u))
	                    return NC_ERR; 


                      //Read v
                      if (!(varb= dataFile.get_var("v")))
	                    return NC_ERR; 
                      if (!varb->set_cur(0,0,0,0))
	                    return NC_ERR; 
                      if (!varb->get(&vG[0],ntime,s_rho,eta_v,xi_rho))
	                    return NC_ERR; 

                      //Read w
                      if (!(varb= dataFile.get_var("omega")))
	                    return NC_ERR; 
                      if (!varb->set_cur(0,0,0,0))
	                    return NC_ERR; 
                      if (!varb->get(&wG[0],ntime,s_w,eta_rho,xi_rho))
	                    return NC_ERR; 




     i1=0;i2=0;i3=0;i4=0;i5=0;
     for(int nt=0;nt<ntime;nt++)
     {


         for (int j=0;j<eta_rho;j++)
         {
             for (int i=0;i<xi_rho;i++)
             {
                  free_surface[t][j][i]=free_surfaceG[i1];
                  i1++;
             }
         }

         for (int k=0;k<s_rho;k++)
         {
              for (int j=0;j<eta_rho;j++)
              {
                   for (int i=0;i<xi_u;i++)
                   {

                      u[t][k][j][i]=uG[i2];
                      i2++;
                   }
              } 
              for (int j=0;j<eta_v;j++)
              { 
                   for (int i=0;i<xi_rho;i++)
                   {
                      v[t][k][j][i]=vG[i4];
                      i4++;
                   }
   
              } 
         }
         for (int k=0;k<s_w;k++)
         {
              for (int j=0;j<eta_rho;j++)
              {
                   for (int i=0;i<xi_rho;i++)
                   {
                      w[t][k][j][i]=wG[i5];
                      i5++;
                   }
              } 
         }
     t++;
    }



     auxi=auxi+ntime;
     dataFile.close();// close Netcdf file
 }




   delete[] free_surfaceG;
   delete[] uG;
   delete[] vG;
   delete[] wG;


}








int ReadingGridFlow()
{


  sprintf(ncfile, "/data/geo/escola/roms_canaryislands/cb_2009_3km_42/roms_avg.0330.nc");
  NcFile dataFile(ncfile,NcFile::ReadOnly);
  if(!dataFile.is_valid())
  return NC_ERR;


  NcDim *ncdim_xi_rho;
  if (!(ncdim_xi_rho = dataFile.get_dim("xi_rho")))
  return NC_ERR;
  xi_rho = ncdim_xi_rho->size();

  NcDim *ncdim_eta_rho;
  if (!(ncdim_eta_rho = dataFile.get_dim("eta_rho")))
  return NC_ERR;
  eta_rho = ncdim_eta_rho->size();

  NcDim *ncdim_s_rho;
  if (!(ncdim_s_rho = dataFile.get_dim("s_rho")))
  return NC_ERR;
  s_rho = ncdim_s_rho->size();

  NcDim *ncdim_xi_u;
  if (!(ncdim_xi_u = dataFile.get_dim("xi_u")))
    return NC_ERR;
  xi_u = ncdim_xi_u->size();
 
  NcDim *ncdim_eta_v;
  if (!(ncdim_eta_v = dataFile.get_dim("eta_v")))
    return NC_ERR;
  eta_v = ncdim_eta_v->size();

  NcDim *ncdim_s_w;
  if (!(ncdim_s_w = dataFile.get_dim("s_w")))
    return NC_ERR;
  s_w= ncdim_s_w->size();

  NcDim *ncdim_time;
  if (!(ncdim_time = dataFile.get_dim("time")))
    return NC_ERR;
  ntime = ncdim_time->size();

  dataFile.close();// close Netcdf file

  xi=xi_u;
  eta=eta_v;
  s=s_w;




  sprintf(ncfileNetwork, "/data/geo/escola/roms_canaryislands/cb_2009_3km_42/cb_2009_3km_grd_smooth.nc");
  NcFile dataNetwork(ncfileNetwork,NcFile::ReadOnly);
  if(!dataNetwork.is_valid())
  return NC_ERR;




  if (!(ncdim_xi_rho = dataNetwork.get_dim("xi_psi")))
  return NC_ERR;
  xi_psi = ncdim_xi_rho->size();


  if (!(ncdim_eta_rho = dataNetwork.get_dim("eta_psi")))
  return NC_ERR;
  eta_psi = ncdim_eta_rho->size();




  NcVar *varb;
  double valor;

  lon_rho = new double *[eta_rho];
  lat_rho = new double *[eta_rho];
  h = new double *[eta_rho];
  mask_rho = new double *[eta_rho];
  lon_psi = new double *[eta_psi];
  lat_psi = new double *[eta_psi];

   for (int i = 0; i < eta_psi; i++)
   {
     lon_psi[i] = new double [xi_psi];
     lat_psi[i] = new double [xi_psi];
   }
   for (int i = 0; i < eta_rho; i++)
   {
     lon_rho[i] = new double [xi_rho];
     lat_rho[i] = new double [xi_rho];
     h[i] = new double [xi_rho];
     mask_rho[i] = new double [xi_rho];
   }



   for (int i = 0; i < xi_rho; i++)
   {
      for (int j = 0; j < eta_rho; j++)
      {


               //Read lon_rho
               if (!(varb= dataNetwork.get_var("lon_rho")))
	             return NC_ERR; 
               if (!varb->set_cur(j, i))
	             return NC_ERR; 
               if (!varb->get(&valor,1,1))
	             return NC_ERR; 

               lon_rho[j][i]=valor;


               //Read lat_rho
               if (!(varb= dataNetwork.get_var("lat_rho")))
	             return NC_ERR; 
               if (!varb->set_cur(j, i))
	             return NC_ERR; 
               if (!varb->get(&valor,1,1))
	             return NC_ERR; 

               lat_rho[j][i]=valor;



               //Read h
               if (!(varb= dataNetwork.get_var("h")))
	             return NC_ERR; 
               if (!varb->set_cur(j, i))
	             return NC_ERR; 
               if (!varb->get(&valor,1,1))
	             return NC_ERR; 

               h[j][i]=valor;


               //Read Mask_rho
               if (!(varb= dataNetwork.get_var("mask_rho")))
	             return NC_ERR; 
               if (!varb->set_cur(j, i))
	             return NC_ERR; 
               if (!varb->get(&valor,1,1))
	             return NC_ERR; 

               mask_rho[j][i]=valor;


      }
   }


   for (int i = 0; i < xi_psi; i++)
   {
      for (int j = 0; j < eta_psi; j++)
      {


               //Read lon_psi
               if (!(varb= dataNetwork.get_var("lon_psi")))
	             return NC_ERR; 
               if (!varb->set_cur(j, i))
	             return NC_ERR; 
               if (!varb->get(&valor,1,1))
	             return NC_ERR; 

               lon_psi[j][i]=valor;

               //Read lat_psi
               if (!(varb= dataNetwork.get_var("lat_psi")))
	             return NC_ERR; 
               if (!varb->set_cur(j, i))
	             return NC_ERR; 
               if (!varb->get(&valor,1,1))
	             return NC_ERR; 

               lat_psi[j][i]=valor;


      }
  }


  dataNetwork.close();// close Netcdf file
  dataFile.close();// close Netcdf file

}





