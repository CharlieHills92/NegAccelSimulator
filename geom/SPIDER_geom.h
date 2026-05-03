#ifndef _SPIDER_GEOM
#define _SPIDER_GEOM

#include "../param/TEST_parameters.h"

////////////////// Geometry of the electrodes  /////////////////////////

/// SPIDER Accelerator
// P. Veltri, IO,  01/2020.


// PG
// here PG start at z=0. B field MUST be adjusted accordingly.

bool solid0_SPIDER( double x, double y, double z )
{
 double A1=-0.8333;
 double b1=0.012;
 double z1=0.006;
 
 double A2=1;
 double b2=0.001;
 double z2=0.009;



/* this part is for multi-beamlet simulation. commend for single beamlet
 if (abs(x)>Lx05*N_aperture_x){  // the If here is to stop aperture repetition in the geom where in the location there are no real beamlets
   if (z<=z2) return true;
   else return false;}
 else{
 double NNNx=floor((x+Lx05)/(2*Lx05));
 x= x-NNNx*2*Lx05;}


 if (abs(y)>Ly05*N_aperture_y){
   if (z<=z2) return true;
   else return false;}
 else{
 double NNNy=floor((y+Ly05)/(2*Ly05));
 y= y-NNNy*2*Ly05;}
*/



 double r=sqrt(x*x+y*y);
 if(z <= z1 && z>=0.0)   return  ((r-(A1*z+b1))>=1e-6); // (r>=A1*z+b1); 
 else if (  z > z1 && z<=z2) return ((r-(A2*z+b2))>=1e-6);
 else return false;

}






// EG

bool solid1_SPIDER( double x, double y, double z )
{
double z3=0.026;
double z2=0.023;
double z1=0.015; 

double A=0.33333;
double b=-0.0011667;

if (abs(x)>Lx05*N_aperture_x){  // the If here is to stop aperture repetition in the geom where in the location there are no real beamlets
   if (z>=z1 && z<=z3) return true;
   else return false;}
 else{
 double NNNx=floor((x+Lx05)/(2*Lx05));
 x= x-NNNx*2*Lx05;}

 if (abs(y)>Ly05*N_aperture_y){
   if (z>=z1 && z<=z3) return true;
   else return false;}
 else{
 double NNNy=floor((y+Ly05)/(2*Ly05));
 y= y-NNNy*2*Ly05;}


double r=sqrt(x*x+y*y);

if(z <= z2 && z>=z1) return ((r-0.0065)>=1e-6) ;
else if(z>z2 && z<=z3 ) return (r-(A*z+b)>=1e-6);
else return false;

//old geom below
//if(z <= z2 && z>=z1) return (r>=0.0065) ;
//else if(z>z2 && z<=z3 ) return r>=A*z+b;
//else return false;


}    





// GG

bool solid2_SPIDER( double x, double y, double z )
{
double z3=0.0735;
double z2=0.0695;
double z1=0.061; 


if (abs(x)>Lx05*N_aperture_x){  // the If here is to stop aperture repetition in the geom where in the location there are no real beamlets
   if (z>=z1 && z<=z3) return true;
   else return false;}
 else{
 double NNNx=floor((x+Lx05)/(2*Lx05));
 x= x-NNNx*2*Lx05;}


 if (abs(y)>Ly05*N_aperture_y){
   if (z>=z1 && z<=z3) return true;
   else return false;}
 else{
 double NNNy=floor((y+Ly05)/(2*Ly05));
 y= y-NNNy*2*Ly05;}




 double r=sqrt(x*x+y*y);

double A=0.125;
double b=-0.0001875;
if (z >z1 && z<z2) return ((r-0.008)>=1e-6) ;
else if(z>=z2 && z<=z3 ) return (r-(A*z+b)>=1e-6);
else return false;
}


#endif

