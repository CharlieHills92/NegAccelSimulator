#include "THCallback.h"




int check_periodicity( ParticleDataBase3D* _pdb, ParticleBase *particle, ParticleP3D *pcur, ParticleP3D *pend, const vector<double>& _periodicity ) {
	
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
        // Round to nearest rather than letting C++ integer division truncate toward zero.
        // The two wrap counters are packed additively as hor*100 + ver*10000, and truncation
        // mis-decodes every mixed-sign pair: a particle one cell along +x and one along -y
        // carries gen = 100-10000 = -9900, which truncates to ver=0, hor=-99 instead of
        // ver=-1, hor=+1. The cell test below then evaluates beamlet_x-99 = -96, fails the
        // >0 condition, and the particle silently stops wrapping for the rest of its life.
        // Rounding recovers the intended pair for every offset the cell test permits
        // (|hor| <= 2, |ver| <= 10), because no legal |hor| reaches the 50 that would make
        // the packing genuinely ambiguous.
        int current_periodic_ver=(int)lround((double)genORI/ver_addition);
        int current_periodic_hor=(int)lround((double)(genORI-ver_addition*current_periodic_ver)/hor_addition);

        if ( beamlet_x+current_periodic_hor<5 && beamlet_x+current_periodic_hor>0 && 
             beamlet_y+current_periodic_ver<16 && beamlet_y+current_periodic_ver>0 ) {
            // Resolve the x and y wraps together and inject exactly ONE replacement.
            //
            // The x block and the y block used to be two independent `if`s, each killing the
            // parent and calling add_particle() with the parent's full IQ. A step that landed
            // outside BOTH bounds -- a diagonal exit, which the transverse sweep makes common
            // -- therefore fired both and injected TWO full-weight particles for one parent,
            // duplicating charge. Neither copy was at the true diagonal neighbour either:
            // one was shifted in x with y still out of bounds, the other vice versa.
            //
            // Accumulating the two shifts and adding one particle conserves IQ exactly and
            // lands a diagonal escape on the correct diagonal cell in a single step. The
            // generation offsets accumulate the same way, so both cell hops stay recorded.
            const double span_x = _periodicity[0]-_periodicity[1];
            const double span_y = _periodicity[2]-_periodicity[3];
            double newx = (*pend)[1];
            double newy = (*pend)[3];
            int gen_shift = 0;

            // x
            if( (*pend)[1] > _periodicity[0] ) {
                bound=1;
                newx -= span_x;
                gen_shift += hor_addition;
            }
            else if( (*pend)[1] < _periodicity[1] ) {
                bound=-1;
                newx += span_x;
                gen_shift -= hor_addition;
            }
            // y
            if( (*pend)[3] > _periodicity[2] ) {
                bound=2;
                newy -= span_y;
                gen_shift += ver_addition;
            }
            else if( (*pend)[3] < _periodicity[3] ) {
                bound=-2;
                newy += span_y;
                gen_shift -= ver_addition;
            }

            if( gen_shift != 0 ) {
                particle->set_status( PARTICLE_PERIODIC );
                ParticleP3D peri( (*pend)[0], newx,(*pend)[2],newy,(*pend)[4],(*pend)[5],(*pend)[6] );
                _pdb->add_particle( particle->IQ(), particle->q()/CHARGE_E, particle->m()/MASS_U,
                                    genORI+gen_shift, peri);
            }
        }
	}
    
    // double delta = sqrt( ((*pend)[1]-(*pcur)[1])*((*pend)[1]-(*pcur)[1])+((*pend)[3]-(*pcur)[3])*((*pend)[3]-(*pcur)[3])+((*pend)[5]-(*pcur)[5])*((*pend)[5]-(*pcur)[5]) );
    // if (debug) logfile << "DEBUG: bound " << bound << " mass " << particle->m() << endl << flush;
    // if (debug) logfile << "DEBUG: delta " <<  delta << " x " << (*pend)[1] << " vx " << (*pend)[2] << " y " << (*pend)[3] << " vy " << (*pend)[4] << " z " << (*pend)[5] << " vz " << (*pend)[6] << endl << flush;
	return bound;
}