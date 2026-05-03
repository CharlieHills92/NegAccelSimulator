#ifndef _GEOM
#define _GEOM


////////////////// Geometry of the electrodes  /////////////////////////

/// ELISE Accelerator, IPP, 
// P. Veltri, IO,  01/2019.


// PG
   //return(x <= 0.0065  && x>0.0 && y>= -0.9992*x+0.013295  && y> 0.007);
   //return( x <= 0.009 && x>0.0065 && y>=0.3997*x+0.0044021);
bool solid0( double x, double y, double z )
{
 

// with PG starting at z=0 // NOTE B field map should be adjusted consequently!
 double z0=0.0;
 //double A1=-0.8413; 
 double A1=-tan(40.0*M_PI/180.0); 
 double b1=0.0123;
 double z1=0.0063;

 double z2=0.0065;

 double A3=1.0;
 double b3=0.0005;
 double z3=0.009;



 if (abs(x)>Lx05*N_aperture_x){  // the If here is to stop aperture repetition in the geom where in the location there are no real beamlets
   if (z<=z3) return true;
   else return false;}
 else{
 double NNNx=floor((x+Lx05)/(2*Lx05));
 x= x-NNNx*2*Lx05;}


 if (abs(y)>Ly05*N_aperture_y){
   if (z<=z3) return true;
   else return false;}
 else{
 double NNNy=floor((y+Ly05)/(2*Ly05));
 y= y-NNNy*2*Ly05;}





 double r=sqrt(x*x+y*y);
 if(z <= z1 && z>=z0)   return  ((r-(A1*z+b1))>1e-6); // (r>=A1*z+b1); 
 else if (  z > z1 && z<=z2) return ((r-0.007)>1e-6);
 else if (  z > z2 && z<=z3) return ((r-(A3*z+b3))>1e-6);
 else return false;

}






// EG

bool solid1( double x, double y, double z )
{
double z3=0.0264;
double z2=0.0194;
double z1=0.015; 

double A=0.1429;
double b=0.0032286;

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


if(z <= z2 && z>=z1) return (r>=0.006) ;
else if(z>z2 && z<=z3 ) return r>=A*z+b;
else return false;
}    

// GG+REPELLER

bool solid2( double x, double y, double z )
{
double z3=0.0654;
double z2=0.0644;
double z1=0.0414; 
double A=1;
double b=-0.0439;


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


if (z >z1 && z<z2) return (r>=0.0075);
else if(z>z2 && z<=z3 ) return (r>=A*z+b);
else return false;
}


////// End of Geometry 
////////////////////////////////////////////////////////////////////////////////





#endif

