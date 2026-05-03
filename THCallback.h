#ifndef THCALLBACK_H_
#define THCALLBACK_H_

#include "particledatabase.hpp"
#include "cross_sections.h"
#include "funct.h"
#include "random.hpp"
#include "globals.h"
#include "constants.hpp"
#include <cmath>
#include <algorithm>
//~ #include "src/pier/matlabfun.h"

using namespace std;

// Add safety constants
#define MAX_GENERATION_DEPTH 5
#define MAX_PARTICLE_COUNT 10000000
#define MIN_VELOCITY_THRESHOLD 1e3
#define MAX_VELOCITY_SCALE SPEED_C

//~ void get_dens_at_z( double z, double & nH2, double & nH0 ) {

    //~ double arr[] = { -0.05,0.,0.005,0.01,0.0155,0.026,0.03,0.035,0.05,0.055,0.062,0.064,0.074,0.08,0.084,0.1,0.15 };
    //~ vector<double> pos( arr, arr+sizeof(arr)/sizeof(arr[0]) );
    //~ double arr2[] = { 2.50E+19,2.50E+19,2.60E+19,2.65E+19,2.50E+19,1.70E+19,1.50E+19,1.41E+19,1.38E+19,1.34E+19,1.11E+19,1.00E+19,5.00E+18,2.50E+18,2.10E+18,2.00E+18,2.00E+18 };
    //~ vector<double> dens( arr2, arr2+sizeof(arr2)/sizeof(arr2[0]) );

//~ };


int check_periodicity( ParticleDataBase3D* _pdb, ParticleBase *particle, ParticleP3D *pcur, ParticleP3D *pend, vector<double>& _periodicity );

class THCallback_test : public TrajectoryHandlerCallback {
public:

    THCallback_test() {}

    virtual ~THCallback_test() {}

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) {
        ParticleP3D *pcur = (ParticleP3D *)( xcur );
        ParticleP3D *pend = (ParticleP3D *)( xend );
        
        //~ cout << "DEBUG -> v1=" << pcur->speed() << " v2=" << pend->speed() << " x1=" << (*pcur)[5] << " x2=" << (*pend)[5] << endl;
        
        // Kill particles with mass less than 0.5 atomic mass and x-coordinate more than 3 mm.
        if( particle->m() > 1.5e-27 && (*pcur)[5] >= 20.0e-3 ) {
            *pend = *pcur;
            particle->set_status( PARTICLE_STRIP );
        }
    }
};



class THCallback_strip : public TrajectoryHandlerCallback {

private:

    Random* _rng;
    ofstream* _debugstream;
    ParticleDataBase3D* _pdb;
    bool _debugprint;
    double _mass;
    vector<double> pos;
    vector<double> dens;
    vector<double> _periodicity;
    bool _periodic;
    
    int _counter;


public:

    THCallback_strip( bool debugprint, ParticleDataBase3D* pdb, double& mass ) {
        vector<double> period_empty;
        period_empty.clear();
        _periodic=false;
        THCallback_strip( debugprint, pdb, mass, period_empty, "densprofiles/MITICA_dens.txt" );
    }

    THCallback_strip( bool debugprint, ParticleDataBase3D* pdb, double& mass, vector<double>& periodicity ) {
        THCallback_strip( debugprint, pdb, mass, periodicity, "densprofiles/MITICA_dens.txt" );
    }

    THCallback_strip( bool debugprint, ParticleDataBase3D* pdb, double& mass, vector<double>& periodicity, string density_filename ) {
        //Random *rng; // Random number generator for two uniform random positions, and three Gaussian velocity components (to get a Maxwell-Boltzmann distribution for the total velocity).
        // if( ibsimu.get_rng_type() == RNG_SOBOL )  //In addition, two more numbers are needed for the cosine distribution and the angles.
        // 	_rng = new QRandom( 0 );
        // else
        // 	_rng = new MTRandom( 0 );
        _rng = new MTRandom( 1 );
        double qx[1];
        _rng->get( qx );
        //~ cout << "DEBUG: random generator initialized\n";
        // DEBUG
        _debugprint=debugprint;
        _pdb=pdb;
        _mass=mass;
        _periodicity=periodicity;
        _periodic=false;
        if (_periodicity.size()>0) {
            _periodic=true;
            logfile << " THC_STRIPPING: THC PERIODICITY IS ON " << endl;
            logfile << " Periodic boundaries: ";
            for (int qq=0; qq<4; qq++) logfile << _periodicity[qq] << " ";
            logfile << endl << flush;
        }
        else {
            if (debug) {
                logfile << " THC_STRIPPING: THC PERIODICITY IS OFF " << endl;
            }
        }

        _counter=0;
        load_density_profile(density_filename, pos, dens);
        // cout << "Immhererad\n";
        // for(uint kk=0; kk<pos.size(); kk++) cout << pos[kk] << "\t" << dens[kk] << endl;
        // cout << "Immhererad\n";
        if( _debugprint ) {
            _debugstream=new ofstream("THCdebug.txt");
            *_debugstream << "zstart(m)" << " ";
            *_debugstream << "zend(m)" << " ";
            *_debugstream << "energy(eV)" << " ";
            *_debugstream << "sigma(m-2)" << " ";
            *_debugstream << "density(m-3)" << " ";
            *_debugstream << "delta(m)" << " ";
            *_debugstream << "mfp(m)" << " ";
            *_debugstream << "rand" << " ";
            *_debugstream << "prob_coll" << " ";
            *_debugstream << "collided" << "\n" << flush;
        }
    }

    virtual ~THCallback_strip() {
        if( _debugprint ) {
            _debugstream->close();
        }
        
    }

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) {
        ParticleP3D *pcur = (ParticleP3D *)( xcur );
        ParticleP3D *pend = (ParticleP3D *)( xend );
        double vel = pcur->speed();
        // double energy = 0.5*particle->m()*vel*vel/1.6e-19;
		double energy = particle->m()*SPEED_C2*(1./sqrt(1.-(vel*vel)/(SPEED_C2))-1.)/CHARGE_E; // in eV
        // cout << "zLOC = " << (*pend)[5] << endl;
        if ( (*pend)[5]>0.007 && (*pend)[6]>0  ) {
            // Reduce current according to stripping of H- probability
            //~ if( particle->m() > 1.5e-27 && energy >= 10. && (particle->get_status() == PARTICLE_OK) ) {
            if( particle->m() > 1.5e-27*_mass && particle->m() < 1.8e-27*_mass && particle->q()<0 && energy >= 10. ) {
                double single_strip;
                double double_strip;
                double sigma = stripping_cross_at_E( energy, _mass, single_strip, double_strip );
                double zc=(*pcur)[5];
                bool ciaone;
                double target_n=density_at_z( zc, 0.3, pos, dens, ciaone );
                //double nu_scat=sigmav*target_n; //
                
                double delta = sqrt( ((*pend)[1]-(*pcur)[1])*((*pend)[1]-(*pcur)[1])+((*pend)[3]-(*pcur)[3])*((*pend)[3]-(*pcur)[3])+((*pend)[5]-(*pcur)[5])*((*pend)[5]-(*pcur)[5]) );
                //double delta = ((*pend)[5]-(*pcur)[5]);
                
                //double dt=delta/vel;
                //double lambda=0.05;
                double lambda=1./(sigma*target_n);
                
                //~ double Prob_collision=1.-exp(-sigma*target_n*delta);
                double Prob_collision=1.-exp(-delta/lambda);
                //~ double Prob_collision=0.01;
                double rand1[1];
                _rng->get( rand1 );
                
                int collided=0;
                
                _counter++;

                if ( rand1[0] < Prob_collision ) {
                    particle->set_status( PARTICLE_STRIP );
                    collided=1;
                }
                
                if ( _debugprint ) {
                    if (collided==1)
                    {
                        *_debugstream << setw(12) << (*pcur)[5] << " ";
                        *_debugstream << setw(12) << (*pend)[5] << " ";
                        *_debugstream << setw(12) << energy << " ";
                        *_debugstream << setw(12) << sigma << " ";
                        *_debugstream << setw(12) << target_n << " ";
                        *_debugstream << setw(12) << delta << " ";
                        *_debugstream << setw(12) << lambda<< " ";
                        *_debugstream << setw(12) << rand1[0] << " ";
                        *_debugstream << setw(12) << Prob_collision << " ";
                        *_debugstream << setw(12) << collided << "\n" << flush;
                    }
                    
                }
            }

            // PERIODIC DOMAIN
            if (_periodic) {
                cout << " ISPERIODIC!" << endl;
                // periodicity vector contains the dimensions of the two periodicities:
                // _periodicity[0]: positive x
                // _periodicity[1]: negative x
                // _periodicity[2]: positive y
                // _periodicity[3]: negative y
                check_periodicity(_pdb, particle, pcur, pend, _periodicity);
            }
        }

    }
        
};

class THCallback_secondaries : public TrajectoryHandlerCallback {

private:

    Random* _rng;
    Random* _Gng;
    ParticleDataBase3D* _pdb;
    double _mass;
    vector<double> pos;
    vector<double> dens;
    bool _periodic;
    vector<double> _periodicity;

    // Helper function for safe velocity scaling
    bool validateAndScaleVelocity(double vel, double minvel, double originalVel[3], double scaledVel[3]) {
        // Validate input velocity
        if (vel < MIN_VELOCITY_THRESHOLD || !std::isfinite(vel)) {
            if (debug) logfile << "WARNING: Invalid velocity in secondary generation: " << vel << endl;
            return false;
        }
        
        // Calculate safe scale factor - prevent energy increase
        double scale_factor = std::min(1.0, minvel/vel);
        
        // Apply velocity scaling with bounds checking
        for (int i = 0; i < 3; i++) {
            scaledVel[i] = scale_factor * originalVel[i];
            
            // Validate scaled velocity
            if (!std::isfinite(scaledVel[i]) || abs(scaledVel[i]) > MAX_VELOCITY_SCALE) {
                if (debug) logfile << "WARNING: Invalid scaled velocity component " << i << ": " << scaledVel[i] << endl;
                scaledVel[i] = 0.0; // Set to rest if invalid
            }
        }
        
        return true;
    }
    
    // Position validation
    bool validatePosition(ParticleP3D *pend) {
        for (int i = 0; i < 7; i++) {
            if (!std::isfinite((*pend)[i])) {
                logfile << "ERROR: Invalid particle position component " << i << " in secondary generation: " << (*pend)[i] << endl;
                return false;
            }
        }
        return true;
    }

public:

    THCallback_secondaries( ParticleDataBase3D* pdb, double& mass ) {
        vector<double> period_empty;
        period_empty.clear();
        _periodic=false;
        THCallback_secondaries( pdb, mass, period_empty, "densprofiles/MITICA_dens.txt" );
    }

    THCallback_secondaries( ParticleDataBase3D* pdb, double& mass, vector<double>& periodicity ) {
        THCallback_secondaries( pdb, mass, periodicity, "densprofiles/MITICA_dens.txt" );
    }

    THCallback_secondaries( ParticleDataBase3D* pdb, double& mass, vector<double>& periodicity, string density_filename ) {
        _pdb=pdb;
        _mass=mass;
        if (debug) cout << " MASS : " << _mass << endl;
        _rng = new MTRandom( 1 );
        _Gng = new MTRandom( 3 );
        double qx[1];
        double Rqx[3];
        _Gng->set_transformation(0, Gaussian_Transformation());
        _Gng->set_transformation(1, Gaussian_Transformation());
        _Gng->set_transformation(2, Gaussian_Transformation());
        _rng->get( qx );
        _Gng->get( Rqx );
        load_density_profile(density_filename,pos,dens);
        _periodicity=periodicity;
        _periodic=false;
        if (_periodicity.size()>0) {
            _periodic=true;
            logfile << " THC_SECONDARIES: THC PERIODICITY IS ON " << endl;
            logfile << " Periodic boundaries: ";
            for (int qq=0; qq<4; qq++) logfile << _periodicity[qq] << " ";
            logfile << endl << flush;
        }
        else {
            if (debug) {
                logfile << " THC_SECONDARIES: THC PERIODICITY IS OFF " << endl;
            }
        }
    }

    virtual ~THCallback_secondaries() {}

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) {
        ParticleP3D *pcur = (ParticleP3D *)( xcur );
        ParticleP3D *pend = (ParticleP3D *)( xend );
        
        // Input validation
        double vel = pcur->speed();
        if (!std::isfinite(vel)) {
            if (debug) logfile << "WARNING: Skipping secondary generation for invalid velocity: " << vel << endl;
            return;
        }
        
        // Validate particle position
        if (!validatePosition(pend)) {
			logfile << "ERROR: Invalid particle position detected. "
					<< "Mass=" << particle->m()
					<< " Charge=" << particle->q()
					<< " Current=" << particle->IQ()
					<< " status=" << particle->get_status()
					<< " x = " << (*pend)[1]
					<< " vx = " << (*pend)[2]
					<< " y = " << (*pend)[3]
					<< " vy = " << (*pend)[4]
					<< " z = " << (*pend)[5]
					<< " vz = " << (*pend)[6]
					<< " xs = " << (*pcur)[1]
					<< " vxs = " << (*pcur)[2]
					<< " ys = " << (*pcur)[3]
					<< " vys = " << (*pcur)[4]
					<< " zs = " << (*pcur)[5]
					<< " vzs = " << (*pcur)[6]
					<< endl << flush;
            return;
        }
        
        // Limit generation depth to prevent infinite chains
        int generation = particle->gen();
        if (generation > MAX_GENERATION_DEPTH) {
            if (debug) logfile << "WARNING: Maximum generation depth reached for particle, stopping secondaries" << endl;
            return;
        }
        
        // Check particle count to prevent memory overflow
        if (_pdb->size() > MAX_PARTICLE_COUNT) {
            logfile << "WARNING: Maximum particle count reached (" << _pdb->size() << "), stopping secondary generation" << endl;
            return;
        }
        
		double energy = particle->m()*SPEED_C2*(1./sqrt(1.-(vel*vel)/(SPEED_C2))-1.)/CHARGE_E; // in eV

        double minvelH = sqrt(2.*10.*CHARGE_E/particle->m()); // velocity corresponding to 10 eV for H/D
        double minvelH2 = sqrt(2.*10.*CHARGE_E/(2*particle->m())); // velocity corresponding to 10 eV for H2/D2
        double minvele = sqrt(2.*0.1*CHARGE_E/(MASS_E)); // velocity corresponding to 0.1 eV for e

        // Validate minimum velocities
        if (!std::isfinite(minvelH) || !std::isfinite(minvelH2) || !std::isfinite(minvele)) {
            logfile << "ERROR: Invalid minimum velocity calculation" << endl;
            return;
        }
                
        double zc=(*pcur)[5];
        bool ciaone = true;
        double target_n=density_at_z( zc, 0.3, pos, dens, ciaone );
        
        // Validate target density
        if (!std::isfinite(target_n) || target_n < 0) {
            if (debug) logfile << "WARNING: Invalid target density at z=" << zc << ": " << target_n << endl;
            return;
        }
        //double nu_scat=sigmav*target_n; //
        
        double delta = sqrt( ((*pend)[1]-(*pcur)[1])*((*pend)[1]-(*pcur)[1])+((*pend)[3]-(*pcur)[3])*((*pend)[3]-(*pcur)[3])+((*pend)[5]-(*pcur)[5])*((*pend)[5]-(*pcur)[5]) );
        
        // Validate step size
        if (!std::isfinite(delta) || delta < 0) {
            if (debug) logfile << "WARNING: Invalid step size in secondary generation: " << delta << endl;
            return;
        }

		if (delta == 0 ) {
			return;
		}
                    
        // Reduce current according to stripping of H- probability
        if ( (*pend)[5]>0.007 && (*pend)[6]>0  ) {
            if ( energy >= 10. && abs(particle->gen()%100)<5 ) {
                if( particle->m() > 1.5e-27*_mass && particle->q() < 0 ) {
                    double single_strip;
                    double double_strip;
                    double sigma_strip = stripping_cross_at_E( energy, _mass, single_strip, double_strip );
					// single_strip = 10*single_strip;
					// double_strip = 10*double_strip;
					// sigma_strip = single_strip+double_strip;
                    double sigma_ioniz = cs_bkg_ionization( energy, _mass );
                    double sigma = sigma_strip+sigma_ioniz;
                    
                    // Validate cross sections
                    if (!std::isfinite(sigma) || sigma < 0) {
                        if (debug) logfile << "WARNING: Invalid cross section at E=" << energy << ": " << sigma << endl;
                        return;
                    }

                    double Prob_collision=1.-exp(-sigma*target_n*delta);
                    double rand1[1];
                    _rng->get( rand1 );
                    
                    if ( rand1[0] < Prob_collision ) {
                        *pend = *pcur;
                        // Null collision
                        double frac_single=single_strip/sigma;
                        double frac_double=(sigma_strip)/sigma;
                        _rng->get( rand1 );
                        //cout << rand1[0] << "\t" << frac_single << "\n";
                        
                        try {
                            if ( rand1[0] <= frac_single ) {
                                // single stripping --> H- diventa neutro
                                particle->set_status( PARTICLE_STRIP );
                                int genORI=particle->gen();
                                _pdb->add_particle( particle->IQ(), 0., 1., genORI+1, *pend); // H0 veloce
                                _pdb->add_particle( particle->IQ(), -1., 1./1836, genORI+1, *pend); // 1e is emitted. Assume it has the same velocity of the original particle 
                            }
                            else if ( rand1[0] <= frac_double ){
                                // double stripping --> H- diventa H+
                                particle->set_status( PARTICLE_STRIP );
                                int genORI=particle->gen();
                                _pdb->add_particle( -particle->IQ(), 1., 1., genORI+1, *pend); // H+ veloce
                                _pdb->add_particle( 2*particle->IQ(), -1., 1./1836, genORI+1, *pend); // 2e are emitted. Assume they have the same velocity of the original particle 
                                
                            }
                            else {
                                // ionizzazione del gas di background -> alla particella non succede nulla. Si genera un H2+ fermo 
                                // particle->set_status( PARTICLE_OK );
                                int genORI=particle->gen();
                                
                                // Safe velocity scaling for H2+ and electron
                                double originalVel[3] = {(*pend)[2], (*pend)[4], (*pend)[6]};
                                double scaledVelH2[3], scaledVelE[3];
                                
                                if (validateAndScaleVelocity(vel, minvelH2, originalVel, scaledVelH2)) {
                                    ParticleP3D seco1( (*pend)[0], (*pend)[1], scaledVelH2[0], (*pend)[3], scaledVelH2[1], (*pend)[5], scaledVelH2[2] );
                                    _pdb->add_particle( -particle->IQ(), 1., 2., genORI+1, seco1); // H2+
                                }
                                
                                if (validateAndScaleVelocity(vel, minvele, originalVel, scaledVelE)) {
                                    ParticleP3D seco2( (*pend)[0], (*pend)[1], scaledVelE[0], (*pend)[3], scaledVelE[1], (*pend)[5], scaledVelE[2] );
                                    _pdb->add_particle( particle->IQ(), -1., 1./1836, genORI+1, seco2); // e
                                }
                            }
                        } catch (const std::exception& e) {
                            logfile << "ERROR: Failed to add secondary particle in H- collision: " << e.what() << endl;
                            // Continue with original particle
                        }
                    }
                }
                // H0 collisions
                else if( particle->m() > 1.5e-27*_mass  && particle->m() < 2.e-27*_mass  && particle->q() == 0 ) {
                    double sigma_strip = cs_proj_ionization_H0( energy, _mass );
                    double sigma_ioniz = cs_bkg_ionization( energy, _mass );
                    double sigma = sigma_strip+sigma_ioniz;
                    
                    // Validate cross sections
                    if (!std::isfinite(sigma) || sigma < 0) {
                        if (debug) logfile << "WARNING: Invalid H0 cross section: " << sigma << endl;
                        return;
                    }

                    double Prob_collision=1.-exp(-sigma*target_n*delta);
                    double rand1[1];
                    _rng->get( rand1 );
                    if ( rand1[0] < Prob_collision ) {
                        *pend = *pcur;
                        // Null collision
                        double frac_strip=sigma_strip/sigma;
                        _rng->get( rand1 );
                        
                        try {
                            if ( rand1[0] <= frac_strip ) {
                                // stripping --> H0 diventa H+
                                particle->set_status( PARTICLE_STRIP );
                                int genORI=particle->gen();
                                _pdb->add_particle( -particle->IQ(), 1., 1., genORI+1, *pend); // H+ veloce
                                _pdb->add_particle( particle->IQ(), -1., 1./1836, genORI+1, *pend); // 1e is emitted. Assume it has the same velocity of the original particle 
                            }
                            else {
                                // ionizzazione del fondo --> H2 diventa H2+, H0 rimane H0
                                int genORI=particle->gen();
                                
                                double originalVel[3] = {(*pend)[2], (*pend)[4], (*pend)[6]};
                                double scaledVelH2[3], scaledVelE[3];
                                
                                if (validateAndScaleVelocity(vel, minvelH2, originalVel, scaledVelH2)) {
                                    ParticleP3D seco1( (*pend)[0], (*pend)[1], scaledVelH2[0], (*pend)[3], scaledVelH2[1], (*pend)[5], scaledVelH2[2] );
                                    _pdb->add_particle( -particle->IQ(), 1., 2., genORI+1, seco1); // H2+
                                }
                                
                                if (validateAndScaleVelocity(vel, minvele, originalVel, scaledVelE)) {
                                    ParticleP3D seco2( (*pend)[0], (*pend)[1], scaledVelE[0], (*pend)[3], scaledVelE[1], (*pend)[5], scaledVelE[2] );
                                    _pdb->add_particle( particle->IQ(), -1., 1./1836, genORI+1, seco2); // e
                                }
                            }
                        } catch (const std::exception& e) {
                            logfile << "ERROR: Failed to add secondary particle in H0 collision: " << e.what() << endl;
                        }
                    }
                }

                // H+ collisions
                else if( particle->m() > 1.5e-27*_mass && particle->m() < 2.e-27*_mass && particle->q() > 0. ) {
                    double sigma_CX = cs_CX_Hp( energy, _mass );
                    double sigma = sigma_CX;
                    
                    // Validate cross section
                    if (!std::isfinite(sigma) || sigma < 0) {
                        if (debug) logfile << "WARNING: Invalid H+ cross section: " << sigma << endl;
                        return;
                    }

                    double Prob_collision=1.-exp(-sigma*target_n*delta);
                    double rand1[1];
                    _rng->get( rand1 );
                    if ( rand1[0] < Prob_collision ) {
                        // CX con il fondo --> H+ diventa H0, H2 diventa H2+
                        *pend = *pcur;
                        particle->set_status( PARTICLE_STRIP ); // TODO mettere PARTICLE_CX
                        int genORI=particle->gen();
                        
                        try {
                            double originalVel[3] = {(*pend)[2], (*pend)[4], (*pend)[6]};
                            double scaledVelH2[3];
                            
                            if (validateAndScaleVelocity(vel, minvelH2, originalVel, scaledVelH2)) {
                                ParticleP3D seco1( (*pend)[0], (*pend)[1], scaledVelH2[0], (*pend)[3], scaledVelH2[1], (*pend)[5], scaledVelH2[2] );
                                _pdb->add_particle( -particle->IQ(), 1., 2., genORI+1, seco1); // H2+
                            }
                            
                            _pdb->add_particle( particle->IQ(), 0., 1., genORI+1, *pend); // H0 veloce
                        } catch (const std::exception& e) {
                            logfile << "ERROR: Failed to add secondary particle in H+ collision: " << e.what() << endl;
                        }
                    }
                }

                // H2+ collisions
                else if( particle->m() > 3.e-27*_mass && particle->m() < 4.e-27*_mass && particle->q() > 0. ) {
                    // TODO - add implementation with same safety measures
                }

                // H0 collisions (molecular)
                else if( particle->m() > 3.e-27*_mass && particle->m() < 4.e-27*_mass && particle->q() == 0. ) {
                    // TODO - add implementation with same safety measures
                }

                // e collisions
                else if( particle->m() < 1e-30 && particle->q() < 0. ) {
                    // TODO - add implementation with same safety measures
                }

            }
        
            // PERIODIC DOMAIN
            if (_periodic) {
                // periodicity vector contains the dimensions of the two periodicities:
                // _periodicity[0]: positive x
                // _periodicity[1]: negative x
                // _periodicity[2]: positive y
                // _periodicity[3]: negative y
                check_periodicity(_pdb, particle, pcur, pend, _periodicity);
            }
        }
    }

};



#endif /* THCALLBACK_H_ */
