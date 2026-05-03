#ifndef _MITICA_GEOM
#define _MITICA_GEOM

//#include "../param/TEST_parameters.h"
#include "../globals.h"


////////////////// Geometry of the electrodes  /////////////////////////

/// MITICA Accelerator, RFX 
// P. Veltri, IO,  01/2019.


// PG
   //return(x <= 0.0065  && x>0.0 && y>= -0.9992*x+0.013295  && y> 0.007);
   //return( x <= 0.009 && x>0.0065 && y>=0.3997*x+0.0044021);
bool solid0_MITICA( double x, double y, double z )
{

 double z0=0.000;
 double A1=-0.82857;
 double b1=0.0128;
 double z1=0.007;

 double z2=0.0074;

 double A3=1.25;
 double b3=-0.00225;
 double z3=0.0078;


 double A4=0.41666;
 double b4=0.00425;
 double z4=0.009;


//utilizzo una condizione che è valida solo per il singolo beamlet...
//perché non voglio avere PG piena a monte

/*
	 if (abs(x)>Lx05*N_aperture_x){  // the If here is to stop aperture repetition in the geom where in the location there are no real beamlets
	   if (z<=z4) return true;
	   else return false;}
	 else{
	 double NNNx=floor((x+Lx05)/(2*Lx05));
	 x= x-NNNx*2*Lx05;}
	
	
	 if (abs(y)>Ly05*N_aperture_y){
	   if (z<=z4) return true;
	   else return false;}
	 else{
	 double NNNy=floor((y+Ly05)/(2*Ly05));
	 y= y-NNNy*2*Ly05;}
*/
	 if (abs(x)>mesh_x){  // the If here is to stop aperture repetition in the geom where in the location there are no real beamlets
	   if (z<=z4) return true;
	   else return false;}
	 else{
	 double NNNx=floor((x+Lx05)/(2*Lx05));
	 x= x-NNNx*2*Lx05;}
	
	
	 if (abs(y)>mesh_y){
	   if (z<=z4) return true;
	   else return false;}
	 else{
	 double NNNy=floor((y+Ly05)/(2*Ly05));
	 y= y-NNNy*2*Ly05;}



 double r=sqrt(x*x+y*y);
 if(z <= z1 && z>=z0)   return  ((r-(A1*z+b1))>1e-6); // (r>=A1*z+b1); 
 else if (  z > z1 && z<=z2) return (r>=0.007);
 else if (  z > z2 && z<=z3) return ((r-(A3*z+b3))>1e-6);
 else if (  z > z3 && z<=z4) return ((r-(A4*z+b4))>1e-6);
 else return false;

}



// EG

bool solid1_MITICA( double x, double y, double z )
{
double z3=0.032;
double z2=0.021;
double z1=0.015; 

double A=0.181818182;//0.16667;
double b=0.002681818;//0.00316667;

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


if(z <= z2 && z>=z1) return (r>=0.0065) ;
else if(z>z2 && z<=z3 ) return r>=A*z+b;
else return false;
}

// AG1

bool solid2_MITICA( double x, double y, double z )
{

double z1=0.120; 
double z2=0.137;


 if (abs(x)>Lx05*N_aperture_x){  // the If here is to stop aperture repetition in the geom where in the location there are no real beamlets
   if (z>=z1 && z<=z2) return true;
   else return false;}
 else{
 double NNNx=floor((x+Lx05)/(2*Lx05));
 x= x-NNNx*2*Lx05;}

 if (abs(y)>Ly05*N_aperture_y){
   if (z>=z1 && z<=z2) return true;
   else return false;}
 else{
 double NNNy=floor((y+Ly05)/(2*Ly05));
 y= y-NNNy*2*Ly05;}

double r=sqrt(x*x+y*y);


if (z >z1 && z<z2) return (r>=0.007);
else return false;
}

// AG2
bool solid3_MITICA( double x, double y, double z )
{

double z1=0.225;
double z2=0.242; 


 if (abs(x)>Lx05*N_aperture_x){  // the If here is to stop aperture repetition in the geom where in the location there are no real beamlets
   if (z>=z1 && z<=z2) return true;
   else return false;}
 else{
 double NNNx=floor((x+Lx05)/(2*Lx05));
 x= x-NNNx*2*Lx05;}

 if (abs(y)>Ly05*N_aperture_y){
   if (z>=z1 && z<=z2) return true;
   else return false;}
 else{
 double NNNy=floor((y+Ly05)/(2*Ly05));
 y= y-NNNy*2*Ly05;}

double r=sqrt(x*x+y*y);


if (z >z1 && z<z2) return (r>=0.007);
else return false;
}


// AG3
bool solid4_MITICA( double x, double y, double z )
{

double z1=0.330;
double z2=0.347; 


 if (abs(x)>Lx05*N_aperture_x){  // the If here is to stop aperture repetition in the geom where in the location there are no real beamlets
   if (z>=z1 && z<=z2) return true;
   else return false;}
 else{
 double NNNx=floor((x+Lx05)/(2*Lx05));
 x= x-NNNx*2*Lx05;}

 if (abs(y)>Ly05*N_aperture_y){
   if (z>=z1 && z<=z2) return true;
   else return false;}
 else{
 double NNNy=floor((y+Ly05)/(2*Ly05));
 y= y-NNNy*2*Ly05;}

double r=sqrt(x*x+y*y);


if (z >z1 && z<z2) return (r>=0.008);
else return false;
}


// AG4
bool solid5_MITICA( double x, double y, double z )
{

double z1=0.435;
double z2=0.452; 


 if (abs(x)>Lx05*N_aperture_x){  // the If here is to stop aperture repetition in the geom where in the location there are no real beamlets
   if (z>=z1 && z<=z2) return true;
   else return false;}
 else{
 double NNNx=floor((x+Lx05)/(2*Lx05));
 x= x-NNNx*2*Lx05;}

 if (abs(y)>Ly05*N_aperture_y){
   if (z>=z1 && z<=z2) return true;
   else return false;}
 else{
 double NNNy=floor((y+Ly05)/(2*Ly05));
 y= y-NNNy*2*Ly05;}

double r=sqrt(x*x+y*y);


if (z >z1 && z<z2) return (r>=0.008);
else return false;
}


// GG

bool solid6_MITICA( double x, double y, double z )
{

double z1=0.540;
double z2=0.557; 


 if (abs(x)>Lx05*N_aperture_x){  // the If here is to stop aperture repetition in the geom where in the location there are no real beamlets
   if (z>=z1 && z<=z2) return true;
   else return false;}
 else{
 double NNNx=floor((x+Lx05)/(2*Lx05));
 x= x-NNNx*2*Lx05;}

 if (abs(y)>Ly05*N_aperture_y){
   if (z>=z1 && z<=z2) return true;
   else return false;}
 else{
 double NNNy=floor((y+Ly05)/(2*Ly05));
 y= y-NNNy*2*Ly05;}

double r=sqrt(x*x+y*y);


if (z >z1 && z<z2) return (r>=0.008);
else return false;
}




////// End of Geometry 
////////////////////////////////////////////////////////////////////////////////

#endif


