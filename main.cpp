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
#include <memory>
#include <streambuf>
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <sys/stat.h>

#include "ManageSimulation_New.h"
#include "error.hpp"
#include "ibsimu.hpp"
//#include "particleiterator_MOD.hpp"
#include "globals.h"
//~ #include "sammy/sammy.h"

using namespace std;

namespace {

bool fileExists(const string& path) {
    ifstream file(path.c_str());
    return file.good();
}

string stemFromPath(const string& path) {
	const size_t slash = path.find_last_of("/\\");
	const size_t start = (slash == string::npos) ? 0 : slash + 1;
	const size_t dot = path.find_last_of('.');

	if (dot != string::npos && dot > start) {
		return path.substr(start, dot - start);
	}

	return path.substr(start);
}

string folderFromPath(const string& path) {
	const size_t slash = path.find_last_of("/\\");
	if (slash == string::npos) {
		return ".";
	}

	if (slash == 0) {
		return "/";
	}

	return path.substr(0, slash);
}

string resolveConfigFilePath(const string& input) {
    if (input.size() >= 5 && input.substr(input.size() - 5) == ".json") {
        return input;
    }

    if (fileExists(input)) {
        return input;
    }

    return input + ".json";
}

string toLowerCopy(const string& value) {
	string normalized = value;
	std::transform(normalized.begin(), normalized.end(), normalized.begin(),
				   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
	return normalized;
}

bool isAbsolutePath(const string& path) {
	return !path.empty() && path[0] == '/';
}

void ensureDirectoryExists(const string& path) {
	if (path.empty()) {
		return;
	}

	string current;
	for (size_t index = 0; index < path.size(); ++index) {
		current += path[index];
		if (path[index] != '/' && index + 1 != path.size()) {
			continue;
		}

		string directory = current;
		while (!directory.empty() && directory[directory.size() - 1] == '/') {
			directory.erase(directory.size() - 1);
		}
		if (directory.empty()) {
			continue;
		}

		if (mkdir(directory.c_str(), 0777) == -1 && errno != EEXIST) {
			throw Error(ERROR_LOCATION, "Could not create directory " + directory + ": " + strerror(errno));
		}
	}
}

string resolveLogFilePath(const string& foldname, const string& configured_log_file, const string& case_tag) {
	string filename = configured_log_file.empty() ? case_tag + "_log.log" : configured_log_file;
	if (!isAbsolutePath(filename)) {
		filename = foldname + "/" + filename;
	}

	const size_t slash = filename.find_last_of('/');
	if (slash != string::npos) {
		ensureDirectoryExists(filename.substr(0, slash));
	}
	return filename;
}

void configureConsoleVerbosity(const string& configured_level) {
	const string level = toLowerCopy(configured_level);

	int verbose_level = 1;
	int warning_level = 1;
	int error_level = 1;
	int debug_general_level = 0;
	int debug_dxf_level = 0;

	if (level == "trace") {
		verbose_level = 2;
		debug_general_level = 2;
		debug_dxf_level = 2;
	} else if (level == "debug") {
		verbose_level = 2;
		debug_general_level = 1;
	} else if (level == "info") {
		verbose_level = 1;
	} else if (level == "warning") {
		verbose_level = 0;
	} else if (level == "error" || level == "critical") {
		verbose_level = 0;
		warning_level = 0;
	}

	ibsimu.set_message_threshold(MSG_VERBOSE, verbose_level);
	ibsimu.set_message_threshold(MSG_WARNING, warning_level);
	ibsimu.set_message_threshold(MSG_ERROR, error_level);
	ibsimu.set_message_threshold(MSG_DEBUG_GENERAL, debug_general_level);
	ibsimu.set_message_threshold(MSG_DEBUG_DXF, debug_dxf_level);
}

bool levelEnablesDebug(const string& configured_level) {
	const string level = toLowerCopy(configured_level);
	return level == "debug" || level == "trace";
}

class TeeStreambuf : public std::streambuf {
public:
	TeeStreambuf(std::streambuf* primary, std::streambuf* secondary)
		: primary_(primary), secondary_(secondary) {}

protected:
	virtual int overflow(int ch) {
		if (ch == traits_type::eof()) {
			return traits_type::not_eof(ch);
		}

		const int primary_result = primary_ ? primary_->sputc(static_cast<char>(ch)) : ch;
		const int secondary_result = secondary_ ? secondary_->sputc(static_cast<char>(ch)) : ch;
		if (primary_result == traits_type::eof() || secondary_result == traits_type::eof()) {
			return traits_type::eof();
		}
		return ch;
	}

	virtual int sync() {
		const int primary_sync = primary_ ? primary_->pubsync() : 0;
		const int secondary_sync = secondary_ ? secondary_->pubsync() : 0;
		return (primary_sync == 0 && secondary_sync == 0) ? 0 : -1;
	}

private:
	std::streambuf* primary_;
	std::streambuf* secondary_;
};

class ScopedStreamRedirect {
public:
	ScopedStreamRedirect(std::ostream& stream, std::streambuf* secondary)
		: stream_(stream), original_(stream.rdbuf()), tee_(original_, secondary) {
		stream_.rdbuf(&tee_);
	}

	~ScopedStreamRedirect() {
		stream_.rdbuf(original_);
	}

private:
	std::ostream& stream_;
	std::streambuf* original_;
	TeeStreambuf tee_;
};

} // namespace

int main( int argc, char **argv )
{
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " <case.json|case_stem> [load_existing]" << endl;
		return 1;
	}

	const string config_file = resolveConfigFilePath(argv[1]);
	if (!fileExists(config_file)) {
		cerr << "Resolved JSON case file does not exist: " << config_file << endl;
		return 1;
	}

	SimulationParameters startup_parameters;
	try {
		startup_parameters.readParametersFromFile(config_file);
	} catch (Error& e) {
		e.print_error_message(cerr, false);
		return 1;
	}

	const string case_tag = stemFromPath(config_file);
	const string foldname = folderFromPath(config_file);
	const string logfilename = resolveLogFilePath(
		foldname,
		startup_parameters.getOutputLoggingStructuredLogFile(),
		case_tag);

	logfile.open( logfilename.c_str() );
	if (!logfile.is_open()) {
		cerr << "Could not open log file: " << logfilename << endl;
		return 1;
	}

	std::unique_ptr<ScopedStreamRedirect> cout_redirect;
	std::unique_ptr<ScopedStreamRedirect> cerr_redirect;
	if (startup_parameters.getOutputLoggingCaptureStdout()) {
		cout_redirect.reset(new ScopedStreamRedirect(cout, logfile.rdbuf()));
		cerr_redirect.reset(new ScopedStreamRedirect(cerr, logfile.rdbuf()));
	}

	configureConsoleVerbosity(startup_parameters.getOutputLoggingConsoleLevel());
	debug = startup_parameters.getOutputLoggingWriteDebugArtifacts() ||
		levelEnablesDebug(startup_parameters.getOutputLoggingConsoleLevel()) ||
		levelEnablesDebug(startup_parameters.getOutputLoggingFileLevel());
	ibsimu.set_message_output(cout);
// logfile << "\tDisplacement grid created. Xdispl " << xdisplacementSG << " Ydispl " << ydisplacementSG << endl; 
	
	logfile << " CONFIG FILE: " << config_file << "\n";
	logfile << " CASE TAG: " << case_tag << "\n";
	logfile << " LOG FILE: " << logfilename << "\n";
	logfile << " OUTPUT FLAGS: summary=" << startup_parameters.getOutputSummaryEnabled()
	       << " plots=" << startup_parameters.getOutputPlotsEnabled()
	       << " data=" << startup_parameters.getOutputDataEnabled()
	       << " vtk=" << startup_parameters.getOutputVTKEnabled() << "\n";

	const uint nscans = 1;
	for ( size_t scan_index=0; scan_index<nscans; scan_index++ ) {

		logfile << "*** Starting simulation " << scan_index+1 << "/" << nscans << " ***\n";
		logfile << "*** Input file: " << config_file << " ***\n" << flush;

		ManageSimulation Simulation( case_tag, foldname );

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
			if (!startup_parameters.getOutputDataEnabled()) {
				cerr << "Cannot use load_existing when outputs.data.enabled is false" << endl;
				logfile << "Cannot use load_existing when outputs.data.enabled is false\n";
				cout_redirect.reset();
				cerr_redirect.reset();
				logfile.close();
				return 1;
			}
		}

		if (dosim) {
				
			double tol=Simulation.get_tolerance();
			uint JEXTiterations=0;

			bool checkEGcurrent = false;
				
			double newJION = Simulation.get_J_ION();
			double targetJ = newJION;

			do {
				JEXTiterations++;
				double z_start = Simulation.get_domain_z_start();
				Simulation.create_geometry(z_start, Simulation.get_domain_z_end(), 1);
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

			logfile << "*** Getting details of particles at end location... " << endl << flush;
			Simulation.particles_end_location();
			logfile << "Done! ***\n" << flush;
			if (startup_parameters.getOutputSummaryEnabled() ||
			    (startup_parameters.getOutputPlotsEnabled() && Simulation.get_stripping()>1)) {
				Simulation.fill_particle_dbs();
			}
			
			const double domain_exit_plane = Simulation.get_domain_z_end() - Simulation.get_MESH_SIZE();
			double zlocsummary = startup_parameters.hasDiagnosticSummaryZPosition()
								 ? startup_parameters.getDiagnosticSummaryZPosition()
								 : domain_exit_plane;
			double emitter_export_z = startup_parameters.hasDiagnosticEmitterExportZPosition()
									  ? startup_parameters.getDiagnosticEmitterExportZPosition()
									  : zlocsummary;

			if (startup_parameters.getOutputSummaryEnabled()) {
				logfile << "*** Start the analysis... " << flush;
				Simulation.analysis(zlocsummary);
			}

			if (startup_parameters.getOutputPlotsEnabled()) {
				Simulation.plot_simulation(argc, argv);
				if (Simulation.get_stripping()>1 && startup_parameters.getDiagnosticWritePerSpeciesPlots()) {
					for (size_t species_index = 0; species_index < particle_kind_count(); ++species_index) {
						Simulation.plot_simulation(argc, argv, int2kind(static_cast<int>(species_index)));
					}
				}
			}
			if (startup_parameters.getOutputSummaryEnabled() || startup_parameters.getOutputPlotsEnabled()) {
				logfile << "Done! ***\n" << flush;
			}

			if (startup_parameters.getOutputSummaryEnabled()) {
				const string finaloutput = Simulation.save_emitter("final_map_outside.txt", emitter_export_z);
				if (!finaloutput.empty()) {
					logfile << "*** Output location of particles saved to file " << finaloutput << " ***\n" << flush;
				}
			}
			logfile << "*** SIMULATION " << scan_index+1 << "/" << nscans <<" COMPLETED ***\n" << flush;

			Simulation.save_results_to_vtk();

		}
		else {

			Simulation.create_geometry(Simulation.get_domain_z_start(), Simulation.get_domain_z_end(), 1);
			logfile << "----- Loading simulation of the full domain -----\n" << flush;
			logfile << "*** Geometry defined ***\n" << flush;
			if (!Simulation.load_simulation()) {
				cerr << "Failed to load saved simulation data" << endl;
				cout_redirect.reset();
				cerr_redirect.reset();
				logfile.close();
				return 1;
			}
			logfile << "*** Simulation loaded ***\n" << flush;
			
			Simulation.particles_end_location();
			if (startup_parameters.getOutputSummaryEnabled() ||
			    (startup_parameters.getOutputPlotsEnabled() && Simulation.get_stripping()>1)) {
				Simulation.fill_particle_dbs();
			}
			
			const double domain_exit_plane = Simulation.get_domain_z_end() - Simulation.get_MESH_SIZE();
			double zlocsummary = startup_parameters.hasDiagnosticSummaryZPosition()
			                     ? startup_parameters.getDiagnosticSummaryZPosition()
			                     : domain_exit_plane;

			if (startup_parameters.getOutputSummaryEnabled()) {
				Simulation.analysis(zlocsummary);
			}
			if (startup_parameters.getOutputPlotsEnabled()) {
				Simulation.plot_simulation(argc, argv);
				if (Simulation.get_stripping()>1 && startup_parameters.getDiagnosticWritePerSpeciesPlots()) {
					for (size_t species_index = 0; species_index < particle_kind_count(); ++species_index) {
						Simulation.plot_simulation(argc, argv, int2kind(static_cast<int>(species_index)));
					}
				}
			}
			if (startup_parameters.getOutputSummaryEnabled() || startup_parameters.getOutputPlotsEnabled()) {
				logfile << "Done! ***\n" << flush;
			}
			logfile << "*** Getting details of particles at end location... " << endl << flush;
			logfile << "Done! ***\n" << flush;
			logfile << "*** SIMULATION " << scan_index+1 << "/" << nscans <<" COMPLETED ***\n" << flush;
		}


	}

	logfile << "*** ALL SIMULATIONS COMPLETED ***";
	cout_redirect.reset();
	cerr_redirect.reset();
	logfile.close();

	

	return 0;

}









