#include <fstream>
#include <iomanip>
#include <limits>

#include <timer.hpp>
#include "random.hpp"
#include "dxf_solid.hpp"
#include "stl_solid.hpp"
#include "stlfile.hpp"
#include "mydxffile.hpp"
#include "trajectorydiagnostics.hpp"
#include "transformation.hpp"
#include "convergence.hpp"
#include "epot_bicgstabsolver.hpp"
#include "particledatabase.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "particlediagplotter.hpp"
#include "geomplotter.hpp"
#include "config.h"
#include "epot_field.hpp"
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <particles.hpp>
#include "./my_diagnostics.h"
#include <sys/stat.h> 
#include <sys/types.h> 
#include <readascii.hpp>

//#include "parameters_BUG.h"
//#include "BUG_geom.h" // this geometry should be fine now!


/*
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include "./my_diagnostics.h"
*/

using namespace std;





double CalcolaAngolo (std::vector<double> &v, std::vector<double> &vz,std::vector<double> &curr)
  {
   double argomento_angolo = 0.0;
   double angolo = 0.0;
   double sumcurr = sumcurrent(curr);

   for( uint32_t i = 0; i < v.size(); i++ ) {
       argomento_angolo += atan2(v[i],vz[i])*curr[i];
       }
   angolo = argomento_angolo/sumcurr;
  return (angolo);
  }

double CalcolaDivergenza (std::vector<double> &v,std::vector<double> &vz,std::vector<double> &curr, double angolo)
{
    double sumcurr = 0.0;
   for( uint32_t i = 0; i < curr.size(); i++ ) {
       sumcurr += curr[i]; 
    }
   double argomento_div = 0.0;

   for( uint32_t i = 0; i < vz.size(); i++ ) {
       double argomento_thetaj2 = atan2(v[i],vz[i]) - angolo;
       double thetaj2 = argomento_thetaj2 * argomento_thetaj2;
       double pesoj = curr[i]/sumcurr;
       // cout << i <<" "<<v[i] <<" "<< vz[i] << " " << angolo << endl;
       if(std::isnan(thetaj2) or std::isnan(pesoj)) {
         thetaj2=0; pesoj=1.0e-20; // if a nan, then put zero and give low weight
         cout << "NaN in a particle n. "<< i << " in divergence calculation!" << endl;
       }
       argomento_div += thetaj2 * pesoj;
       thetaj2 = 0.0;
       pesoj = 0.0; 
    }
   double divergenza = sqrt(argomento_div);

 return (divergenza);
}


double CalcolaDivergenza2 (std::vector<double> &v,std::vector<double> &vz,std::vector<double> &curr, double angolo)
{
uint32_t n=v.size();
std::vector<double> theta;
double divergence;
  for( uint32_t i = 0; i < n; i++ ) {
    theta.push_back(atan2(v[i],vz[i]));
  }
// Not weighted avergae (on current!)
divergence=stdev(theta);

return divergence;
}







double sumcurrent (std::vector<double> &curr)
{
 double sumcurr = 0.0;
   for( uint32_t i = 0; i < curr.size(); i++ ) sumcurr += curr[i]; 
 return (sumcurr);
}



double stdev (std::vector<double> &div_vect){
   double mean=0.0;
   double sd=0.0;

   for( int  n = 0; n < int(div_vect.size()); n++ ){
       mean+=div_vect[n];
    }
   mean=mean/div_vect.size();

   double var = 0;
   for( int n = 0; n < int(div_vect.size()); n++ )
   {
      var += (div_vect[n] - mean) * (div_vect[n]- mean);
   }
   var /= div_vect.size();
   sd = sqrt(var);
return sd;
}



double Average(std::vector<double> &v,std::vector<double> &wei)
  {
   double ave = 0.0;
   
   double nan_count = 0.0;
   
   double wtot=0.0;
   
   for( uint32_t i = 0; i < v.size(); i++ ) {

       if(std::isnan(v[i])){
         cout << "NaN in element n. "<< i << " in average calculation!" << endl;
         nan_count++;}
       else{
         ave += v[i]*wei[i];
         wtot +=wei[i];
       };
    }
   
  ave = ave/wtot;
  return (ave);
  }




double Clearance(Geometry &geom,std::vector<double> &x,std::vector<double> &y,double zloc,double Lx05)
  { // TODO to be done for multi beamlet simulation and for asymmetric or no centered aperturef this assumes round apertures, centered in (0,0). Alternative: divide by beamlet, then find x0, y0 and radius of each beamlet and evaluate the distance from that radius, after ubtraction of x0,y0.
  //double xmax=0.0; double xmin=0.0; //find x0, y0 and radius of aperture 
  
double r=Lx05;
double xmax;
  double res=0.00002; // m, resolution used to find the radius 
  for( uint32_t i = 0; i < 10000; i++ ){ // looks for and edge moving along x
       xmax=x[0]+double(i*res);
      //cout <<"x: "<< xmax << "bool: " << solid1(xmax,y[0],zloc) << endl;
      //if(solid1(xmax,y[0],zloc)){  // withou standar Ibsimu function
        if(geom.inside(Vec3D(xmax,y[0],zloc))>0){ // edge found
         r=sqrt(xmax*xmax+y[0]*y[0]);
         //cout<< "edge found"<<endl;
         break;
       }
  }

//now calculate the clearance particle after particle
     double rpart;
     double clea=999.0;
   for( uint32_t i = 0; i < x.size(); i++ ) {
     rpart=sqrt(x[i]*x[i]+y[i]*y[i]);
     clea=min(clea,r-rpart);
   }
   clea=max(0.0,clea);
  //cout <<"z,r,clear: "<< zloc << " "<< r <<" "<< clea<< endl;
  return(clea);
   }


double Clearance_electrode(Geometry &geom, ParticleDataBase3D &pdb,uint32_t nsolid,std::vector<trajectory_diagnostic_e> & diagnostics,TrajectoryDiagnosticData &tdata, double Lx05,double Ly05){
  std::vector<double> x[1];
  std::vector<double> y[1];
  double z;
  double zmin=999.0;
  double zmax=0.0;
  double res=0.0001; // m, resolution used to find the radius 
  for( uint32_t i = 0; i < 200000; i++ ){ // looks for the z position of the solid
  z=i*res;
  if(geom.inside(Vec3D(Lx05,Ly05,z))==nsolid){
  zmin=min(z,zmin);
  break;}
  }

  for( uint32_t i = 0; i < 200000; i++ ){ // looks for the z position of the solid
  z=zmin+i*res;
  if(geom.inside(Vec3D(Lx05,Ly05,z))!=nsolid)  
  break;
  zmax=max(z,zmax);
  }
//

// TODO tomorrow load all particles (real point traked) in a volume V zmin<V<zmax and check the poitn one by one with geom.inside(). instead of calling many trajectories_at_plane

  // now calculate clearance inside the solid
  res=0.00025; 
  uint32_t nstep=uint32_t(round((zmax-zmin)/res))+1;
  double clea=Lx05;
  double min_clea=clea;
  
  for( uint32_t i = 0; i < nstep; i++ ){ // looks for and edge moving along x
   // save x and y at diagnostic plane for clearance
   x[0].clear();
   y[0].clear();
   pdb.trajectories_at_plane( tdata, AXIS_Z,zmin+i*res, diagnostics );
   if(tdata.traj_size()>0){
      for( uint32_t i = 0; i < tdata.traj_size(); i++) {
        x[0].push_back(tdata(i,0));
        y[0].push_back(tdata(i,1));
      }
       clea=Clearance(geom,x[0],y[0],zmin+i*res,Lx05);
       min_clea=min(clea,min_clea);
   }
  //cout <<"z: "<<zmin+i*res<<" clea: " <<clea << " min clea: "<< min_clea<< endl;
  }
  return(min_clea); // minimum clearance inside a solid
   }



double max_vec(std::vector<double> &v){
  double mymax = v[0]; //let, first element is the smallest one
  for(uint32_t i = 1; i < v.size(); i++){  //start iterating from the second element
    if(v[i] > mymax){
       mymax = v[i];
    }
  }
return(mymax);
}

double min_vec(std::vector<double> &v){
  double mymin = v[0]; //let, first element is the smallest one
  for(uint32_t i = 1; i < v.size(); i++){  //start iterating from the second element
    if(v[i] < mymin){
       mymin = v[i];
    }
  }
return(mymin);
}




void savePot(EpotField &epot, Geometry &geom, string &output_name)
  {
  cout <<" Writing Potential on txt file...  \n";

   string filename = output_name + "_pot3D.txt";
   ofstream pout( filename.c_str() );
double xloc, yloc, zloc;
for(uint32_t i = 0; i < geom.size(0); i++){
xloc=geom.origo(0)+geom.h()*i;
  for(uint32_t j = 0; j < geom.size(1); j++){
  yloc=geom.origo(1)+geom.h()*j;
    for(uint32_t k = 0; k < geom.size(2); k++){
    zloc=geom.origo(2)+geom.h()*k;
      Vec3D pos3d(xloc,yloc,zloc);
      pout << xloc*1000<< " "<< yloc*1000 << " "<< zloc*1000 << " "<< epot(pos3d) << " "<< "\n";
      //cout << xloc*1000<< " "<< yloc*1000 << " "<< zloc*1000 << " "<< epot(pos3d) << " "<< "\n";
    }
  }
}
  cout <<" OK!  \n";
}



void saveEfield(EpotEfield &Efield, Geometry &geom, string &output_name)
  {
  cout <<" Writing Potential on txt file...  \n";

   string filename =  output_name + "_Efield3D.txt";
   ofstream pout( filename.c_str() );
double xloc, yloc, zloc;
for(uint32_t i = 0; i < geom.size(0); i++){
xloc=geom.origo(0)+geom.h()*i;
  for(uint32_t j = 0; j < geom.size(1); j++){
  yloc=geom.origo(1)+geom.h()*j;
    for(uint32_t k = 0; k < geom.size(2); k++){
    zloc=geom.origo(2)+geom.h()*k;
      Vec3D pos3d(xloc,yloc,zloc);
      pout << xloc*1000<< " "<< yloc*1000 << " "<< zloc*1000 << " "<< Efield(pos3d) << " "<< "\n";
      //cout << xloc*1000<< " "<< yloc*1000 << " "<< zloc*1000 << " "<< epot(pos3d) << " "<< "\n";
    }
  }
}
  cout <<" OK!  \n";
}




void save2Dmap(EpotField &epot,MeshScalarField &scharge ,Geometry &geom, string &output_name){

cout <<" Writing 2D slice of Potential and rho on txt file...  \n";

   string filename = output_name + "_pot2D.txt";
   ofstream pout( filename.c_str() );
double xlocal, ylocal, zlocal;
double resolutio=0.5;
uint32_t zlim;
zlim=uint32_t(0.01/geom.h())+1; // # nodes in z, only export up to z=8mm

for(uint32_t i = 0; i < geom.size(0)/resolutio; i++){
  xlocal=geom.origo(0)+geom.h()*i*resolutio;
  ylocal=geom.origo(1)+geom.size(1)*geom.h()/2.0; // cut at half y plane
    for(uint32_t k = 0; k < zlim/resolutio; k++){
    zlocal=geom.origo(2)+geom.h()*k*resolutio;
      Vec3D pos3d(xlocal,ylocal,zlocal);
      pout << xlocal*1000<< " "<< ylocal*1000 << " "<< zlocal*1000 << " "<< epot(pos3d) << " "<< scharge(pos3d) << "\n";
//      cout << zlim <<xlocal*1000<< " "<< ylocal*1000 << " "<< zlocal*1000 << " "<< epot(pos3d) << " "<< "\n";
    }
}
  cout <<" OK!  \n";

  pout.close();
}
















/////

void SavePart(ParticleDataBase3D &pdb, double z_diagn,string &output_name,TrajectoryDiagnosticData &tdata,std::vector<trajectory_diagnostic_e> & diagnostics){

// Save data of all particles at diagnostic plane.
cout << "Saving ALL Particles data at Diagnostic plane z= "<< z_diagn << "m ..."<<endl;

std::vector<double> no_part[1];

pdb.trajectories_at_plane( tdata, AXIS_Z,z_diagn,diagnostics);

//string filename = "./"+ output_name + "_" + "Particles.txt";

string filename =  output_name;

ofstream pout( filename.c_str() );

pout << "# x y z vx vy vz cur mass" << "\n"<<std::flush ;

	if(tdata.traj_size()>0){
   		for( uint32_t i = 0; i < tdata.traj_size(); i++) {
		//uint32_t l = no_part[0][i];
		//cout << "l "<<tdata(l,0) << " i "<< tdata(i,0) <<endl; ma perche usava l invece di i ????
    		pout << tdata(i,0) << " "
		<< tdata(i,1) << " "
		<< tdata(i,2) << " "
		<< tdata(i,3) << " "
		<< tdata(i,4) << " "
		<< tdata(i,5) << " "
		<< tdata(i,6) << " "
		<< tdata(i,7) << " "
		<< "\n"<<std::flush ;
		}

	} // end if
	
	pout.close();
	cout <<" OK!  " << endl;
} // end funct SavePart

///////





//save trajectories (TODO, number of particles printed is not always correct)
void saveTraj(ParticleDataBase3D &pdb, string &output_name,uint32_t div){
string traj_name;
string dir_name ="./"+ output_name + "/Traj/";

ParticleP3D part;

uint32_t N=round(pdb.size()/div);


  if (N>pdb.size()){
      cout << "Error: try to print " << N << "particle, while database only contains:" << pdb.size() << endl;
      }
//int div=floor(pdb.size()/N);


// Creating a directory 
    if (mkdir(dir_name.c_str(), 0777) == -1) 
        cerr << "Error when creating directory for Trajectories:  " << strerror(errno) << endl; 
    else
        cout << "Directory for Trajectory saving created" << endl; 

 // T Kalvas method for creating folder (example in web site: "Thirteenth example: Rough/Fine mesh simulation"
    // Make directory
    //std::string cmd = "mkdir " + outdir;
    //FILE *fp = popen( cmd.c_str(), "r" );
    //pclose( fp );
    //outdir += "/";




// Now prints
  cout << "Printing Trajectories of "<< N <<" particles on file (ASCII) in the folder: " << dir_name << endl; 
  
	if(N>0){ 
  		for(uint32_t i = 0; i < N; i++){
		traj_name=dir_name+"/t_" + std::to_string(i) + ".txt";
		ofstream pout( traj_name.c_str() );
			for(uint32_t j =0; j<pdb.traj_size(i*div); j++){
			Particle3D &pp=pdb.particle(i);
			part=pdb.trajectory_point(i*div,j);
			pout << part.location()[0] <<" "<<part.location()[1] << " "
			<<part.location()[2] << " "<<part.velocity()[0] << " "
			<<part.velocity()[1] << " "<<part.velocity()[2] << " "
			<< pp.m()/MASS_U<< " " <<pp.IQ() << "\n";
    			}
  		}
  	}	

}

//////////////////////



/////////////////////
void loadEmitter(ParticleDataBase3D *pdb, std::string &emit_filename,double sign_v, double charge_factor, double mass_color){
// read emitter from file

// if needed reverse velocities (sign_v)
//se necessario attenua la carica
// or change mass by 1/10000 only to change colors to particles in the plots. mass color should be 0 < m_color< 100 to avoid affecting the mass too much (>1%)

// Input particles
cout << "Loading the emitter from external file: " << emit_filename<< "\n";
ReadAscii din( emit_filename, 8 );
cout << "Reading " << din.rows() << " particles\n";
// Go through all read particles
// x y z vx vy vz cur mass
for( size_t i = 0; i < din.rows(); i++ ) {
    
    //double m_file  = din[6][i];
    //double t  = din[2][i];
    double x_file  = din[0][i];
    double vx_file = din[3][i];
    double y_file  = din[1][i];
    double z_file  = din[2][i];
    double vy_file = din[4][i];
    double vz_file = din[5][i];
    double I_file  = din[6][i];
    double m_file = din[7][i];
    //cout <<"n. "<<i << " x: "<<x_file << " y: "<<y_file <<" vx: "<<vx_file <<" vy: "<<vy_file <<" vz: "<<vz_file << "m: "<< m_file << endl;
    pdb->add_particle( I_file*charge_factor, -1.0, m_file+m_file*mass_color/10000, ParticleP3D(0.0,x_file,sign_v*vx_file,y_file,sign_v*vy_file,z_file,sign_v*vz_file) );
}
cout << "Emitter Loading: OK!\n";

}






//void plotTraj()



