
#include "funct.h"



particle_kind identify_particle_species( double mass, double charge, double ION_MASS ) {

	particle_kind kind=PARTICLE_WRONG;

	if ( mass<1e-30 ) {
		kind=PARTICLE_E; // electrons
	}
	else if ( (mass>ION_MASS*1.5e-27) & (mass<ION_MASS*1.8e-27) & (charge<0) ) {
		kind=PARTICLE_HM; // neg ions
	}
	else if ( (mass>ION_MASS*1.5e-27) & (mass<ION_MASS*1.8e-27) & (charge>0) ) {
		kind=PARTICLE_HP; // H+
	}
	else if ( (mass>ION_MASS*1.5e-27) & (mass<ION_MASS*1.8e-27) & (charge==0) ) {
		kind=PARTICLE_H0; // H0
	}
	else if ( (mass>2*ION_MASS*1.5e-27) & (mass<2*ION_MASS*1.8e-27) & (charge>0) ) {
		kind=PARTICLE_H2P; // H2+
	}
	else if ( (mass>2*ION_MASS*1.5e-27) & (mass<2*ION_MASS*1.8e-27) & (charge>0) ) {
		kind=PARTICLE_H20; // H20
	}
	else {
		kind=PARTICLE_WRONG; // wrong def
	}

	return kind;

}

string get_particle_name(particle_kind pk) {
	string name;
	switch (pk) {
		case PARTICLE_HM:
			name="HM";
			break;
		case PARTICLE_H0:
			name="H0";
			break;
		case PARTICLE_HP:
			name="HP";
			break;
		case PARTICLE_H2P:
			name="H2P";
			break;
		case PARTICLE_H20:
			name="H20";
			break;
		case PARTICLE_E:
			name="E";
			break;
		case PARTICLE_WRONG:
			name="WRONG";
			break;
		case PARTICLE_ALL:
			name="ALL";
			break;
		default:
			name="";
			break;
	}
	return name;
}

int get_particle_int(particle_kind pk) {
	int num;
	switch (pk) {
		case PARTICLE_HM:
			num=0;
			break;
		case PARTICLE_H0:
			num=1;
			break;
		case PARTICLE_HP:
			num=2;
			break;
		case PARTICLE_H2P:
			num=3;
			break;
		case PARTICLE_H20:
			num=4;
			break;
		case PARTICLE_E:
			num=5;
			break;
		case PARTICLE_WRONG:
			num=-1;
			break;
		case PARTICLE_ALL:
			num=-100;
			break;
		default:
			num=-100;
			break;
	}
	return num;
}


particle_kind int2kind(int num) {
	particle_kind pk;
	switch (num) {
		case 0:
			pk=PARTICLE_HM;
			break;
		case 1:
			pk=PARTICLE_H0;
			break;
		case 2:
			pk=PARTICLE_HP;
			break;
		case 3:
			pk=PARTICLE_H2P;
			break;
		case 4:
			pk=PARTICLE_H20;
			break;
		case 5:
			pk=PARTICLE_E;
			break;
		case -1:
			pk=PARTICLE_WRONG;
			break;
		case -100:
			pk=PARTICLE_ALL;
			break;
		default:
			pk=PARTICLE_ALL;
			break;
	}

	return pk;
}

double density_at_z( double zc, double pressure, vector<double>& pos, vector<double>& dens, bool & ciaone ) {
	
	double target_n=_interp1(pos, dens, zc, ciaone);
	
	return target_n;
	
}


void load_density_profile( string filename, vector<double>& pos, vector<double>& dens ) {
	ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Error opening file!" << std::endl;
    }

    // pos = new vector<double>;
	// dens = new vector<double>;
	
    string line;

    while (getline(inputFile, line)) {
        stringstream ss(line);
        double value1, value2;
        char delimiter;

        // Try comma-separated first
        if (ss >> value1 >> delimiter >> value2 && delimiter == ',') {
            pos.push_back(value1);
            dens.push_back(value2);
        } else {
            // Try tab-separated or space-separated
            stringstream ss2(line);
            if (ss2 >> value1 >> value2) {
                pos.push_back(value1);
                dens.push_back(value2);
            } else {
                cerr << "Error parsing line: " << line << endl;
            }
        }
    }
//	for(uint kk=0; kk<pos->size(); kk++) cout << (*pos)[kk] << "\t" << (*dens)[kk] << endl;

    inputFile.close();
}

double interp1(double x1,double x2, double y1, double y2, double xi){
	return y1+(xi-x1)*(y2-y1)/(x2-x1);
}

double _interp1(vector<double> x, vector<double> y, double xi, bool & found){
	found=false;
	double yi=0;
	if(xi==x[x.size()-1]){
			   found=true; yi=y[x.size()-1];std::cout << yi <<" fatto!\n";
	}
    
	size_t j=1;
    if(xi>=x[0] && xi<=x[x.size()-1]) {
	  while((found==false) && (j<x.size())) {
	    if(xi<x[j]){
		  found=true;
		  yi=interp1(x[j-1],x[j],y[j-1],y[j],xi);
	    }
        j++;
	  }
	}
	if (!found) {
      std::cout << "Warning: value out of range for 1D interpolation!\n";
      std::cout <<"x: "<<x[x.size()-1] << ", xi: "<<xi;
      std::cout <<", x-xi: "<<x[x.size()-1]-xi <<". \n";
      //exit(1);
	}
    //  std::cout << "max x: "<< x[x.size()-1]*1000 <<" xi "<< xi*1000 <<" yi "<< yi*1000 <<" !\n";
    return yi;
}

