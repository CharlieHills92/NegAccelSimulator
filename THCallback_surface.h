////////////////// NOT USED. STILL TO BE DEBUGGED 20240307 //////////////////////



#ifndef THCALLBACK_SURFACE_H_
#define THCALLBACK_SURFACE_H_

#include "particledatabase.hpp"
#include <vector>
#include "globals.h"

using namespace std;

class THCallback_surface : public TrajectoryEndCallback {
private:
    Geometry &_geom;


public:

    THCallback_surface( Geometry &geom )
        : _geom(geom) {
        size_t nsolids = _geom.number_of_solids();
        if (debug) logfile << "DEBUG: Initializing : " << nsolids << endl << flush;
        vector<vector<size_t>>* colls = new vector<vector<size_t>>;
        
        for ( size_t ss=0; ss<nsolids+7; ss++ ) {
            vector<size_t> vpart;
            colls->push_back(vpart);
        }
    }

    virtual ~THCallback_surface() {}

    virtual void operator()( ParticleBase *particle, class ParticleDataBase *pdb ) {

        Particle3D *p3d = (Particle3D *)( particle );
        Vec3D loc = p3d->location();
        Vec3D vel = p3d->velocity();
        double E = 0.5*p3d->m()*vel.ssqr()/pp->q();
        double P = p3d->IQ()*E;



       // Make secondaries based on location and energy
        if( fabs(loc[0]) < 10.1e-3 && loc[1] > 0.0 && 
            loc[2] > 34.9e-3 && loc[2] < 55.1e-3 && E > 1e3 ) {

            // Get normal
            Vec3D normal = _geom.surface_normal( loc );

            // Adjust location off the surface
            loc += 0.01*_geom.h()*normal;

            // Randomize velocity and direction
            double x[3];
            _rand.get( x );
            double mass = 1.0/1836.00;
            double speed = sqrt( 2.0*x[0]*CHARGE_E/(mass*MASS_U) );

            // Find tangents
            Vec3D tang1 = normal.arb_perpendicular();
            Vec3D tang2 = cross( normal, tang1 );
            tang1.normalize();
            tang2.normalize();

            // Build velp in natural coordinates
            double azm_angle = 2.0*M_PI*x[1];
            double pol_angle = asin( sqrt(x[2]) );
            Vec3D velp( speed*cos(pol_angle),  speed*sin(pol_angle), 0.0 );
            Transformation t;
            t.rotate_x( azm_angle );
            velp = t.transform_vector( velp );

            // Convert to surface coordinates
            Vec3D vel = velp[0]*normal + velp[1]*tang1 + velp[2]*tang2;

            ParticleDataBase3D *pdb3d = (ParticleDataBase3D *)( pdb );
            pdb3d->add_particle( p3d->IQ(), -1.0, mass, ParticleP3D( 0.0, 
                                                                     loc[0], vel[0], 
                                                                     loc[1], vel[1], 
                                                                     loc[2], vel[2] ) );
        }
    }

};


#endif
