#include "THCallback.h"




int check_periodicity( ParticleDataBase3D* _pdb, ParticleBase *particle, ParticleP3D *pcur, ParticleP3D *pend, vector<double>& _periodicity ) {
	
	int bound=0;

    // cout << "z position: " << (*pend)[5] << "\n";
    // int x;
    // cin >>x;
	// consider only particles after the PG

	if( (*pend)[5]>0.009 && particle->get_status()==PARTICLE_OK ) {


        int beamlet_x=3;
        int beamlet_y=5;

        int ver_addition=10000;
        int hor_addition=100;

        int genORI=particle->gen();
        int current_periodic_ver=(int)genORI/ver_addition;
        int current_periodic_hor=(int)(genORI-ver_addition*current_periodic_ver)/hor_addition;

        if ( beamlet_x+current_periodic_hor<5 && beamlet_x+current_periodic_hor>0 && 
             beamlet_y+current_periodic_ver<16 && beamlet_y+current_periodic_ver>0 ) {
            // positive x
            if( (*pend)[1] > _periodicity[0] ) {
                bound=1;
                particle->set_status( PARTICLE_PERIODIC );
                double newx=(*pend)[1]-(_periodicity[0]-_periodicity[1]);
                ParticleP3D peri( (*pend)[0], newx,(*pend)[2],(*pend)[3],(*pend)[4],(*pend)[5],(*pend)[6] );
                _pdb->add_particle( particle->IQ(), particle->q()/CHARGE_E, particle->m()/MASS_U, genORI+hor_addition, peri);
            }
            // negative x
            else if( (*pend)[1] < _periodicity[1] ) {
                bound=-1;
                particle->set_status( PARTICLE_PERIODIC );
                double newx=(*pend)[1]+(_periodicity[0]-_periodicity[1]);
                ParticleP3D peri( (*pend)[0], newx,(*pend)[2],(*pend)[3],(*pend)[4],(*pend)[5],(*pend)[6] );
                _pdb->add_particle( particle->IQ(), particle->q()/CHARGE_E, particle->m()/MASS_U, genORI-hor_addition, peri);
            }
            // positive y
            if( (*pend)[3] > _periodicity[2] ) {
                bound=2;
                particle->set_status( PARTICLE_PERIODIC );
                int genORI=particle->gen();
                double newy=(*pend)[3]-(_periodicity[2]-_periodicity[3]);
                ParticleP3D peri( (*pend)[0], (*pend)[1],(*pend)[2],newy,(*pend)[4],(*pend)[5],(*pend)[6] );
                _pdb->add_particle( particle->IQ(), particle->q()/CHARGE_E, particle->m()/MASS_U, genORI+ver_addition, peri);
            }
            // negative y
            else if( (*pend)[3] < _periodicity[3] ) {
                bound=-2;
                particle->set_status( PARTICLE_PERIODIC );
                int genORI=particle->gen();
                double newy=(*pend)[3]+(_periodicity[2]-_periodicity[3]);
                ParticleP3D peri( (*pend)[0], (*pend)[1],(*pend)[2],newy,(*pend)[4],(*pend)[5],(*pend)[6] );
                _pdb->add_particle( particle->IQ(), particle->q()/CHARGE_E, particle->m()/MASS_U, genORI-ver_addition, peri);
            }
        }
	}
    
    // double delta = sqrt( ((*pend)[1]-(*pcur)[1])*((*pend)[1]-(*pcur)[1])+((*pend)[3]-(*pcur)[3])*((*pend)[3]-(*pcur)[3])+((*pend)[5]-(*pcur)[5])*((*pend)[5]-(*pcur)[5]) );
    // if (debug) logfile << "DEBUG: bound " << bound << " mass " << particle->m() << endl << flush;
    // if (debug) logfile << "DEBUG: delta " <<  delta << " x " << (*pend)[1] << " vx " << (*pend)[2] << " y " << (*pend)[3] << " vy " << (*pend)[4] << " z " << (*pend)[5] << " vz " << (*pend)[6] << endl << flush;
	return bound;
}