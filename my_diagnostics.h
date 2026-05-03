#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include "ibsimu.hpp"
#include "geometry.hpp"
#include "epot_field.hpp"
#include "epot_efield.hpp"
#include "meshscalarfield.hpp"


#ifndef _MY_DIAGNOSTICS
#define _MY_DIAGNOSTICS


double CalcolaAngolo (std::vector<double> &,std::vector<double> &,std::vector<double> &);

double CalcolaDivergenza (std::vector<double> &,std::vector<double> &,std::vector<double> &, double);

double CalcolaDivergenza2 (std::vector<double> &,std::vector<double> &,std::vector<double> &, double);

double sumcurrent (std::vector<double> &);

double stdev (std::vector<double> &);

double Average (std::vector<double> &,std::vector<double> &); // weighted average

double max_vec(std::vector<double> &);

double min_vec(std::vector<double> &);

/// check Celarance with all solid at a given z.
double Clearance(Geometry &,std::vector<double> &, std::vector<double> &, double, double); 

/// check clearance inside a specific solid.
double Clearance_electrode(Geometry &,ParticleDataBase3D &,uint32_t, std::vector<trajectory_diagnostic_e> &,TrajectoryDiagnosticData &, double, double); 


/// Save the potential as txt file x(mm),y(mm),z(mm),Pot(V)
void savePot(EpotField &, Geometry &, std::string &);

/// Save the E field map as txt file x(mm),y(mm),z(mm),Ex,Ey,Ez
void saveEfield(EpotEfield &, Geometry &, std::string &);

/// Save a 2D map xz of pot and rho in meniscus area: x,y,z,pot,rho
void save2Dmap(EpotField &,MeshScalarField & ,Geometry &, std::string &);


/// Save all particle on a given plane
void SavePart(ParticleDataBase3D &, double, std::string &,TrajectoryDiagnosticData &,std::vector<trajectory_diagnostic_e> &);


/// Save tracks on a file
void saveTraj(ParticleDataBase3D &, std::string &, uint32_t);



//load emitter form file
void loadEmitter(ParticleDataBase3D *, std::string &, double,double,double);

#endif
