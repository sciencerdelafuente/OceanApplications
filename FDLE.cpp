/**************************************************************************************/  
//  Computes the FDLE (https://link.aps.org/doi/10.1103/PhysRevE.104.065111)
//  Postprocessing program. Input file: Output data from OceanFDLE.cpp:
//  initial and final particle positions at different depths 
//  given in (xi,eta) coordinates (Arakawa C-grid)
/**************************************************************************************/  
//  Code written by Rebeca de la Fuente :  2020                                 
/**************************************************************************************/  
#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string> 

# define PI  3.141592653589793238462643383279502884L 

using namespace std;


int main(void)
{

  

  double dato;
  ifstream archivo;  // objeto de la clase ofstream
  archivo.open("InputFile");



  double Delta=0.5;

  int num,numx,numy;
  double auxiliarx,auxiliary;
  auxiliarx=(331-2)/Delta;
  auxiliary=(533-2)/Delta;
  numx=auxiliarx;
  numy=auxiliary;
  num=auxiliarx*auxiliary;


   int **contador,**domain;
   double **xi,**eta, **xi0,**eta0,**time,**lon0,**lat0;
   contador = new int* [numx];
   xi = new double* [numx];
   eta = new double* [numx];
   xi0 = new double* [numx];
   eta0 = new double* [numx];
   time = new double* [numx];
   domain = new int* [numx];
   lon0 = new double* [numx];
   lat0 = new double* [numx];

   for(int j=0; j<numx; j++)
   {
      contador[j] = new int [numy];
      xi[j] = new double [numy];
      eta[j] = new double [numy];
      xi0[j] = new double [numy];
      eta0[j] = new double [numy];
      time[j] = new double [numy];
      domain[j] = new int [numy];
      lon0[j] = new double [numy];
      lat0[j] = new double [numy];
   } 




  
  int *indicex,*indicey;
  indicex = new int [num];
  indicey = new int [num];





  int cont=0;
  double ix=2;
  for(int i=0;i<numx;i++)
  {
  double iy=2;
      for(int j=0;j<numy;j++)
      {
          contador[i][j]=cont;
          indicex[cont]=i;
          indicey[cont]=j;
      iy=iy+Delta;
      cont++;
      }
  ix=ix+Delta;
  }



  for(int i=0;i<numx;i++)
  {
      for(int j=0;j<numy;j++)
      {
          domain[i][j]=0;
      }
  }



    int ContadorDatos;
    double x0,y0,xt,yt,timet,loni,lati;

   
    while(!archivo.eof())
    {

         archivo >> dato;
         ContadorDatos=dato;
         archivo >> dato;
         x0=dato;
         archivo >> dato;
         y0=dato;
         archivo >> dato;
         archivo >> dato;
         loni=dato;
         archivo >> dato;
         lati=dato;
         archivo >> dato;
         xt=dato;
         archivo >> dato;
         yt=dato;
         archivo >> dato;
         archivo >> dato;
         archivo >> dato;
         archivo >> dato;
         archivo >> dato;
         timet=dato;

         if(!archivo.eof())
         {

         lon0[indicex[ContadorDatos]][indicey[ContadorDatos]]=loni;
         lat0[indicex[ContadorDatos]][indicey[ContadorDatos]]=lati;
         xi0[indicex[ContadorDatos]][indicey[ContadorDatos]]=x0;
         eta0[indicex[ContadorDatos]][indicey[ContadorDatos]]=y0;
         xi[indicex[ContadorDatos]][indicey[ContadorDatos]]=xt;
         eta[indicex[ContadorDatos]][indicey[ContadorDatos]]=yt;
         time[indicex[ContadorDatos]][indicey[ContadorDatos]]=timet;
         domain[indicex[ContadorDatos]][indicey[ContadorDatos]]=1;

         }

    }



   double L1,L2,det,traza,llM;
   double RR[2][2],a[2][2];
   double LyapunovTime,Lyapunovd,Lyapunov;
   int aux;
   int numf=0;




         for(int i=1;i<numx-1;i++)
         {
             for(int j=1;j<numy-1;j++)
             {

                 aux=domain[i][j]*domain[i-1][j]*domain[i][j-1]*domain[i+1][j]*domain[i][j+1];
                 if(aux!=0)
                 {

                  RR[0][0]=((xi[i+1][j]-xi[i-1][j])/(Delta));
                  RR[1][0]=((eta[i+1][j]-eta[i-1][j])/(Delta));

                  RR[0][1]=((xi[i][j+1]-xi[i][j-1])/(Delta));
                  RR[1][1]=((eta[i][j+1]-eta[i][j-1])/(Delta));


                  a[0][0]=RR[0][0]*RR[0][0]+RR[1][0]*RR[1][0];
                  a[0][1]=RR[0][0]*RR[0][1]+RR[1][0]*RR[1][1];
                  a[1][0]=a[0][1];
                  a[1][1]=RR[1][1]*RR[1][1]+RR[0][1]*RR[0][1];


                  traza=a[0][0]+a[1][1];
                  det=a[0][0]*a[1][1]-a[1][0]*a[0][1];
                  L1=(traza+sqrt(traza*traza-4*det))*0.5;
                  L2=(traza-sqrt(traza*traza-4*det))*0.5;
                  L1=fabs(L1);
                  L2=fabs(L2);
                  llM=L1;
                  if(L2>L1){llM=L2;}
                  llM=sqrt(llM);

                  LyapunovTime=log(llM)/(time[i][j]);
                   if(time[i][j]>0 && time[i-1][j]>0 && time[i+1][j]>0 && time[i][j-1]>0 && time[i][j+1]>0)
                   {
                     cout<<xi0[i][j]<<" "<<eta0[i][j]<<" "<<LyapunovTime<<" "<<log(llM)<<" "<<time[i][j]<<endl;
                     numf++;
                   }

                  }
             }
          }














}









