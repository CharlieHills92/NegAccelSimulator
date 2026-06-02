/*
 * main.cpp
 *
 *  Created on: Mar 19, 2020
 *      Author: carlo
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/stat.h> 
#include <sys/types.h> 

#include "ManageSimulation_New.h"
#include "ScanManager.h"
#include "readinputs.h"
//#include "particleiterator_MOD.hpp"
#include "globals.h"
//~ #include "sammy/sammy.h"

using namespace std;

int main( int argc, char **argv )
{

	char* input_cmd = argv[1];
	string scan_file_tag;
	stringstream ss;
	ss << input_cmd;
	ss >> scan_file_tag;

	//string scan_file_name = scan_file_tag+".scn";
	string foldname=scan_file_tag;
	
    if (mkdir(foldname.c_str(), 0777) == -1) cerr << "Warning :  " << strerror(errno) << endl; 
    else cout << "Directory created"; 

	string logfilename=foldname+"/"+scan_file_tag+"_log.log";

	logfile.open( logfilename.c_str() );
// logfile << "\tDisplacement grid created. Xdispl " << xdisplacementSG << " Ydispl " << ydisplacementSG << endl; 
	
	logfile << " SCAN TAG: " << scan_file_tag << "\n";

	// Try to use new ScanManager first, fallback to old system if needed
	vector<string> input_file_tags;
	
	// Check if .scn file exists to determine which system to use
	string scn_file = scan_file_tag + ".scn";
	ifstream scn_check(scn_file);
	bool use_new_scan = scn_check.good();
	scn_check.close();
	
	if (use_new_scan) {
		logfile << "Using ScanManager for file: " << scn_file << "\n";
		ScanManager scanManager;
		
		if (scanManager.loadScanFile(scn_file)) {
			input_file_tags = scanManager.createInputFileFromScan();
			logfile << "Created " << input_file_tags.size() << " input files using ScanManager\n";
		} else {
			logfile << "ScanManager failed, falling back to original method\n";
			input_file_tags = create_input_file_from_scan(scan_file_tag);
		}
	} else {
		logfile << "Using original scan system\n";
		input_file_tags = create_input_file_from_scan(scan_file_tag);
	}
	
	uint nscans = input_file_tags.size();
		
	for ( size_t scan_index=0; scan_index<nscans; scan_index++ ) {

		logfile << "*** Starting simulation " << scan_index+1 << "/" << nscans << " ***\n";
		logfile << "*** Input file: " << input_file_tags[scan_index] << ".inp ***\n" << flush;
		string input_file_tag = input_file_tags[scan_index];

		ManageSimulation Simulation( input_file_tag, foldname );

		bool dosim=true;
		int dosimnum=0;
		if (argc > 2) {
			char* value = argv[2];
			string valstr;
			stringstream ss2;
			ss2 << value;
			ss2 >> valstr;
			dosimnum = atoi( valstr.c_str() );
		}
		if ( dosimnum == 1 ) {
			dosim=false;
	 		ibsimu.message(1) << "	LOADING THE SIMULATION \n";
		}

		if (dosim) {
				
			double tol=Simulation.get_tolerance();
			uint JEXTiterations=0;

			bool checkEGcurrent = false;
				
			double newJION = Simulation.get_J_ION();
			double targetJ = newJION;

			do {
				JEXTiterations++;
				double z_start = 0.;
				if (Simulation.get_accelerator_type() == 1U) z_start = -0.003;
				Simulation.create_geometry(z_start, Simulation.get_domain_z_size(), 1);
				logfile << "----- Simulation of the full domain -----\n" << flush;
				logfile << "*** Geometry defined ***\n" << flush;
				Simulation.add_Bfield();
				logfile << "*** Bfield added ***\n" << flush;
				logfile << "Tolerance "<< tol << "\n" << flush;

				Simulation.set_J_ION(newJION);
				logfile << "*** IT# " << JEXTiterations << " - J_ION_SET = " << newJION << " A/m2 ***\n" << flush;
				Simulation.define_pdb();
				logfile << "*** Emitter defined ***\n" << flush;
				logfile << "*** Launching the simulation... ***\n" << flush;
				Simulation.run_simulation(true);
				logfile << "*** Simulation completed! ***\n" << flush;
			
				double extsimJ = 0.;
				double startsimJ = Simulation.get_J_ION();
				logfile << " Get JSIM\n" << flush;
				checkEGcurrent = Simulation.check_EGext(targetJ, extsimJ);
				logfile << "*** J_ION_SET = " << startsimJ << " A/m2 - J_EXT = "
						<< extsimJ << " A/m2 - J_REQUIRED = " << targetJ << flush;

				if(!checkEGcurrent) {
						newJION=startsimJ/extsimJ*targetJ;
						logfile << "\n --> OUTSIDE " << tol*100. << "% after " << JEXTiterations << " iterations ***\n" << flush;
				}
			}
			while ( (!checkEGcurrent) & (JEXTiterations<5) );
			if (checkEGcurrent) logfile << " --> WITHIN " << tol*100. << "% AFTER " << JEXTiterations << " ITERATIONS ***\n" << flush;
			else logfile << " --> OUTSIDE " << tol*100. << "% AFTER " 
						 << JEXTiterations << " iterations ***\n" << flush;

			logfile << "*** Start the analysis... " << flush;
			logfile << "*** Getting details of particles at end location... " << endl << flush;
			Simulation.particles_end_location();
			logfile << "Done! ***\n" << flush;
			Simulation.fill_particle_dbs();
			
			// Use the actual domain exit plane instead of a hardcoded value.
			// 0.565m is only valid for MITICA/MTF (567mm domain); SPIDER domain is ~85mm.
			double zlocsummary = Simulation.get_domain_z_size() - Simulation.get_MESH_SIZE();

			Simulation.analysis(zlocsummary);
			Simulation.plot_simulation(argc, argv);
			if (Simulation.get_stripping()>1) {
				Simulation.plot_simulation(argc, argv, PARTICLE_HM);
				Simulation.plot_simulation(argc, argv, PARTICLE_H0);
				Simulation.plot_simulation(argc, argv, PARTICLE_HP);
				Simulation.plot_simulation(argc, argv, PARTICLE_H2P);
				Simulation.plot_simulation(argc, argv, PARTICLE_H20);
				Simulation.plot_simulation(argc, argv, PARTICLE_E);
			}
			logfile << "Done! ***\n" << flush;
			
			// Create individual simulation summary for ALL particles with error handling
			try {
				Simulation.create_individual_simulation_summary(scan_index, scan_file_tag, 
					zlocsummary, PARTICLE_ALL);
				logfile << " Individual summary for ALL particles done! ***\n" << flush;
			} catch (const std::exception& e) {
				logfile << " ERROR in individual summary for ALL particles: " << e.what() << "\n" << flush;
			}
			
			// // Add to scan-level beam properties summary for ALL particles with error handling
			// try {
			// 	Simulation.add_to_scan_beam_properties_summary(scan_index, scan_file_tag,
			// 		foldname, scan_file_tag, zlocsummary, PARTICLE_ALL);
			// 	logfile << " Scan beam properties summary for ALL particles done! ***\n" << flush;
			// } catch (const std::exception& e) {
			// 	logfile << " ERROR in scan beam properties summary for ALL particles: " << e.what() << "\n" << flush;
			// }
			
			// // Add to scan-level grid power summary for ALL particles with error handling
			// try {
			// 	Simulation.add_to_scan_grid_power_summary(scan_index, scan_file_tag,
			// 		foldname, scan_file_tag, PARTICLE_ALL);
			// 	logfile << " Scan grid power summary for ALL particles done! ***\n" << flush;
			// } catch (const std::exception& e) {
			// 	logfile << " ERROR in scan grid power summary for ALL particles: " << e.what() << "\n" << flush;
			// }
			
			// if (Simulation.get_stripping()>1) {
			// 	// Create individual summaries for each particle species with error handling
			// 	try {
			// 		Simulation.create_individual_simulation_summary(scan_index, scan_file_tag, 
			// 			zlocsummary, PARTICLE_HM);
			// 		Simulation.create_individual_simulation_summary(scan_index, scan_file_tag, 
			// 			zlocsummary, PARTICLE_H0);
			// 		Simulation.create_individual_simulation_summary(scan_index, scan_file_tag, 
			// 			zlocsummary, PARTICLE_HP);
			// 		Simulation.create_individual_simulation_summary(scan_index, scan_file_tag, 
			// 			zlocsummary, PARTICLE_H2P);
			// 		Simulation.create_individual_simulation_summary(scan_index, scan_file_tag, 
			// 			zlocsummary, PARTICLE_H20);
			// 		Simulation.create_individual_simulation_summary(scan_index, scan_file_tag, 
			// 			zlocsummary, PARTICLE_E);
			// 		logfile << " Individual summaries for all species done! ***\n" << flush;
			// 	} catch (const std::exception& e) {
			// 		logfile << " ERROR in individual summaries for species: " << e.what() << "\n" << flush;
			// 	}
				
			// 	// Add to scan-level beam properties summaries for each species with error handling
			// 	try {
			// 		Simulation.add_to_scan_beam_properties_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, zlocsummary, PARTICLE_HM);
			// 		Simulation.add_to_scan_beam_properties_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, zlocsummary, PARTICLE_H0);
			// 		Simulation.add_to_scan_beam_properties_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, zlocsummary, PARTICLE_HP);
			// 		Simulation.add_to_scan_beam_properties_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, zlocsummary, PARTICLE_H2P);
			// 		Simulation.add_to_scan_beam_properties_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, zlocsummary, PARTICLE_H20);
			// 		Simulation.add_to_scan_beam_properties_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, zlocsummary, PARTICLE_E);
			// 		logfile << " Scan beam properties summaries for all species done! ***\n" << flush;
			// 	} catch (const std::exception& e) {
			// 		logfile << " ERROR in scan beam properties summaries for species: " << e.what() << "\n" << flush;
			// 	}
				
			// 	// Add to scan-level grid power summaries for each species with error handling
			// 	try {
			// 		Simulation.add_to_scan_grid_power_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, PARTICLE_HM);
			// 		Simulation.add_to_scan_grid_power_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, PARTICLE_H0);
			// 		Simulation.add_to_scan_grid_power_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, PARTICLE_HP);
			// 		Simulation.add_to_scan_grid_power_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, PARTICLE_H2P);
			// 		Simulation.add_to_scan_grid_power_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, PARTICLE_H20);
			// 		Simulation.add_to_scan_grid_power_summary(scan_index, scan_file_tag,
			// 			foldname, scan_file_tag, PARTICLE_E);
			// 		logfile << " All particle species summaries done! ***\n" << flush;
			// 	} catch (const std::exception& e) {
			// 		logfile << " ERROR in scan grid power summaries for species: " << e.what() << "\n" << flush;
			// 	}
			// }

			string finaloutput = Simulation.save_emitter("final_map_outside.txt", zlocsummary);
			logfile << "*** Output location of particles saved to file " << finaloutput << " ***\n" << flush;
			logfile << "*** SIMULATION " << scan_index+1 << "/" << nscans <<" COMPLETED ***\n" << flush;

			Simulation.save_results_to_vtk();

		}
		else {

			Simulation.create_geometry(0., Simulation.get_domain_z_size(), 1);
			logfile << "----- Loading simulation of the full domain -----\n" << flush;
			logfile << "*** Geometry defined ***\n" << flush;
			Simulation.load_simulation();
			logfile << "*** Simulation loaded ***\n" << flush;
			
			Simulation.particles_end_location();
			Simulation.fill_particle_dbs();
			
			// Use the actual domain exit plane instead of a hardcoded value.
			double zlocsummary = Simulation.get_domain_z_size() - Simulation.get_MESH_SIZE();

			Simulation.analysis(zlocsummary);
			Simulation.plot_simulation(argc, argv);
			if (Simulation.get_stripping()>1) {
				Simulation.plot_simulation(argc, argv, PARTICLE_HM);
				Simulation.plot_simulation(argc, argv, PARTICLE_H0);
				Simulation.plot_simulation(argc, argv, PARTICLE_HP);
				Simulation.plot_simulation(argc, argv, PARTICLE_H2P);
				Simulation.plot_simulation(argc, argv, PARTICLE_H20);
				Simulation.plot_simulation(argc, argv, PARTICLE_E);
			}
			logfile << "Done! ***\n" << flush;
			logfile << "*** Getting details of particles at end location... " << endl << flush;
			logfile << "Done! ***\n" << flush;
			
			// Create individual simulation summary for ALL particles with error handling
			try {
				Simulation.create_individual_simulation_summary(scan_index, scan_file_tag, 
					zlocsummary, PARTICLE_ALL);
				logfile << " Individual summary for ALL particles done! ***\n" << flush;
			} catch (const std::exception& e) {
				logfile << " ERROR in individual summary for ALL particles: " << e.what() << "\n" << flush;
			}
			logfile << "*** SIMULATION " << scan_index+1 << "/" << nscans <<" COMPLETED ***\n" << flush;
		}


	}

	logfile << "*** ALL SIMULATIONS COMPLETED ***";
	logfile.close();

	

	return 0;

}









