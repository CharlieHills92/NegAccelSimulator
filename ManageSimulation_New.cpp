/*
 * ManageSimulation_New.cpp
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored ManageSimulation)
 */

#include "ManageSimulation_New.h"

#include "globals.h"
#include "ibsimu.hpp"
#include "error.hpp"
#include "epot_bicgstabsolver.hpp"
#include "epot_mgsolver.hpp"
#include "geometry.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "meshscalarfield.hpp"
#include "particledatabase.hpp"
#include "funct.h"
#include "THCallback.h"
#include "THCallback_surf_EAMCC.h"
#include "StrippingUtils.h"

#include <iostream>

using namespace std;

class ForcedPot : public CallbackFunctorB_V {
    public:
        ForcedPot() {}
        ~ForcedPot() {}

        virtual bool operator()( const Vec3D &x ) const {
            // return( x[2] < 0.2e-3 && x[0]*x[0]+x[1]*x[1] > 7e-3*7e-3 );
            return( x[2] < 7e-3 && (abs(x[0])>0.01 || abs(x[1])>0.01) );
            // return( x[2] < 7e-3 );
        }
};

ManageSimulation::ManageSimulation(const std::string& scan_name, const std::string& foldername) {
    try {
        initializeComponents(scan_name, foldername);
        initializeIbsimu();
        
        ibsimu.message(1) << endl << "*** SIMULATION INITIALIZED ***" << endl << endl;
        
    } catch (const exception& e) {
        throw Error(ERROR_LOCATION, "Simulation initialization failed: " + string(e.what()));
    }
}

ManageSimulation::~ManageSimulation() {
    // Cleanup is handled automatically by unique_ptr destructors
}

void ManageSimulation::initializeComponents(const std::string& scan_name, const std::string& foldername) {
    // Initialize managers in dependency order
    parameters = std::unique_ptr<SimulationParameters>(new SimulationParameters());

    const string input_file = foldername + "/" + scan_name + ".json";
    cout << "INPUT FILENAME: " << input_file << endl;

    // Read parameters from the generated JSON case file before configuring output paths.
    parameters->readParametersFromFile(input_file);

    fileManager = std::unique_ptr<FileManager>(new FileManager(
        foldername,
        parameters->getOutputSummaryDirectory(),
        parameters->getOutputPlotsDirectory(),
        parameters->getOutputDataDirectory(),
        parameters->getOutputVTKDirectory()));
    geometryManager = std::unique_ptr<GeometryManager>(new GeometryManager());
    fieldManager = std::unique_ptr<FieldManager>(new FieldManager());
    particleManager = std::unique_ptr<ParticleManager>(new ParticleManager());
    diagnosticsManager = std::unique_ptr<DiagnosticsManager>(new DiagnosticsManager());
    
    // Set file tag after directories are configured.
    fileManager->setFileTag(scan_name);
    diagnosticsManager->setDensityProfileFilename(parameters->getStrippingDensityProfile());
}

void ManageSimulation::initializeIbsimu() {
    // Initialize IBSIMU settings
    ibsimu.set_message_threshold(MSG_VERBOSE, 1);
    const bool monte_carlo_enabled = parameters &&
        (parameters->getIncludeStripping() > 0 || parameters->getIncludeSurfaceCollisions() > 0);
    ibsimu.set_thread_count((debug || monte_carlo_enabled) ? 1 : 8);
}

void ManageSimulation::ResetSimulation() {
    if (geometryManager) geometryManager->resetGeometry();
    if (fieldManager) fieldManager->resetFields();
    if (particleManager) particleManager->resetParticles();
}

void ManageSimulation::create_geometry() {
    geometryManager->createGeometry(*parameters, fileManager->getGeomFile());

    if (parameters->getOutputVTKEnabled() && parameters->getOutputVTKExportGeometry()) {
        const string vtk_base = fileManager->getVTKFolder() + fileManager->getFileTag() + "_geometry";
        export_geometry_to_vtk(vtk_base);
    }
}

void ManageSimulation::create_geometry(double z_start, double z_end, double meshsize_multiplier) {
    geometryManager->createGeometry(*parameters, fileManager->getGeomFile(), 
                                   z_start, z_end, meshsize_multiplier);

    if (parameters->getOutputVTKEnabled() && parameters->getOutputVTKExportGeometry()) {
        const string vtk_base = fileManager->getVTKFolder() + fileManager->getFileTag() + "_geometry";
        export_geometry_to_vtk(vtk_base);
    }
}

void ManageSimulation::add_Bfield() {
    string bfield_folder = FieldManager::getBFieldFolder(*parameters);
    fileManager->setBFieldFolder(bfield_folder);
    fileManager->setBFieldFn(bfield_folder + "/");
    
    fieldManager->addMagneticField(*parameters, bfield_folder);
}

void ManageSimulation::define_pdb() {
    particleManager->defineParticleDatabase(*geometryManager->getGeometry(), *parameters);
}

void ManageSimulation::define_pdb(const std::string& emittername, bool backstream) {
    particleManager->defineParticleDatabase(*geometryManager->getGeometry(), *parameters,
                                           emittername, backstream);
}

void ManageSimulation::run_simulation() {
    run_simulation(false);
}

void ManageSimulation::run_simulation(bool pdbincycle) {
    if (debug) logfile << "DEBUG: Building geometrical surfaces... " << flush;
    
    Geometry* geometry = geometryManager->getGeometry();
    if (!geometry) {
        throw Error(ERROR_LOCATION, "No geometry available for simulation");
    }
    
    if (true) {
        geometry->build_mesh();
        if (debug) logfile << "MESH... " << flush;
    }
    if (true) {
        geometry->build_surface();
        if (debug) logfile << "SURFACES... " << flush;
    }
    if (debug) logfile << "Done!\n" << flush;

    logfile << "Defining solver, plasma, epot, scharge, scharge_ave... " << flush;
    
    // Select solver based on MGSOLVER parameter (0=BiCGSTAB, 1=MG solver)
    EpotSolver* solver = nullptr;
    uint solver_type = parameters->getMGSolver();
    uint meniscus_type = parameters->getShieldModel();

    EpotBiCGSTABSolver* bicgstab_solver = nullptr;
    EpotMGSolver* mg_solver = nullptr;
    
    InitialPlasma initp( AXIS_Z, 7e-3 );
    ForcedPot force;

    if (solver_type == 0) {
        logfile << "Using BiCGSTAB solver... " << flush;
        bicgstab_solver = new EpotBiCGSTABSolver(*geometry);
        if (meniscus_type == 1U) {
            ibsimu.message(1) << " Setting meniscus shield model in BiCGSTAB solver.\n";
            // bicgstab_solver->set_shield_plasma(parameters->getTPositive(), parameters->getUPlasma());
            // bicgstab_solver->set_initial_plasma(parameters->getUPlasma(), &initp);
            bicgstab_solver->set_initial_plasma(0, &initp);
        } else {
            ibsimu.message(1) << " Setting nsimp plasma model in BiCGSTAB solver.\n";
            bicgstab_solver->set_nsimp_initial_plasma( &initp );
        }
        // bicgstab_solver->set_forced_potential_volume( 0.0, &force );
        solver = bicgstab_solver;
    } else {
        logfile << "Using MG solver... " << flush;
        mg_solver = new EpotMGSolver(*geometry);
        if (meniscus_type == 1U) {
            // ibsimu.message(1) << " Setting meniscus shield model in MG solver.\n";
            // mg_solver->set_shield_plasma(parameters->getTPositive(), parameters->getUPlasma());
            // mg_solver->set_initial_plasma(parameters->getUPlasma(), &initp);
            mg_solver->set_initial_plasma(0, &initp);
        } else {
            ibsimu.message(1) << " Setting nsimp plasma model in MG solver.\n";
            mg_solver->set_nsimp_initial_plasma( &initp );
        }
        // mg_solver->set_forced_potential_volume( 0.0, &force );
        solver = mg_solver;
    }
    
    if (debug) logfile << "Done!" << endl << flush;

    // solver->set_gnewton( true ); // enabled by default in BiCGSTAB, not present for MG solver
    
    // Create field objects
    EpotField* epot = new EpotField(*geometry);
    MeshScalarField* scharge = new MeshScalarField(*geometry);
    MeshScalarField scharge_ave(*geometry);
    
    if (debug) logfile << "Done!" << endl << flush;
    
    logfile << "Defining Bfield... " << flush;
    // Get magnetic field from field manager
    MeshVectorField* bfield = fieldManager->getMagnetic();
    
    // Create zero magnetic field if none exists
    MeshVectorField* zero_bfield = nullptr;
    if (bfield == nullptr) {
        // Create zero magnetic field using the mesh and fout parameters
        bool fout[3] = {true, true, true};  // Enable all components
        zero_bfield = new MeshVectorField(*geometry, fout);
        bfield = zero_bfield;  // Use zero field for calculations
    }
    logfile << "Done!\n" << flush;
    
    logfile << "Defining Efield... " << flush;
    EpotEfield* efield = new EpotEfield(*epot);
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield->set_extrapolation(efldextrpl);
    logfile << "Done!\n" << flush;

    logfile << "Defining particle database... " << flush;
    ParticleDataBase3D* pdb = particleManager->getParticles();
    logfile << "Done!\n" << flush;

    // Magnetic field suppression (if available)
    double usup = 10.0; 
    // Note: NPlasmaBfieldSuppression may need to be instantiated if header is found    
    NPlasmaBfieldSuppression psup( *epot, usup );
    pdb->set_bfield_suppression( &psup );
    
    ibsimu.message(1) << "*** SIMULATION READY ***" << endl;

    double rho_h = pdb->get_rhosum(); // moved outside the main loop
    double rho_tot = rho_h;
    cout << " RHO_H = " << rho_h << "\n";

    uint n_iterations = parameters->getIterations();
    
    // Run main simulation loop
    for (uint i = 0; i < n_iterations; ++i) {
        ibsimu.message(1) << "Start iteration # " << i+1 << "/" << n_iterations << "... ";
        logfile << " Start iteration # " << i+1 << "/" << n_iterations << "...\n" << flush;
        
        if (pdbincycle) {
            // Recreate particle database if needed
            particleManager->defineParticleDatabase(*geometry, *parameters);
            pdb = particleManager->getParticles();
        }

        if (i == 1) {
            // Set up plasma conditions for positive ions
            if (meniscus_type == 1U) {
                ibsimu.message(1) << " Setting meniscus shield model in solver.\n";
                solver->set_shield_plasma(parameters->getTPositive(), parameters->getUPlasma());
            } else {
                ibsimu.message(1) << " Setting nsimp plasma model in solver.\n";
                std::vector<double> Ei, rhoi;
                Ei.push_back( parameters->getTPositive() ); // Temperature of positive ions
                rhoi.push_back( 0.5*rho_h );
                double rhop = rho_tot - rho_h*0.5;
                solver->set_nsimp_plasma(rhop, parameters->getUPlasma(), rhoi, Ei);

                std::cout << " rhop = " << rhop << "\n";
                std::cout << " rhoi = ";
                for (size_t idx = 0; idx < rhoi.size(); ++idx) {
                    std::cout << rhoi[idx] << " ";
                }
                std::cout << "\n";
                std::cout << " Uplasma = " << parameters->getUPlasma() << "\n";
                std::cout << " Ei = ";
                for (size_t idx = 0; idx < Ei.size(); ++idx) {
                    std::cout << Ei[idx] << "\n ";
                }
            }
            // std::cout << "\nPress ENTER to continue...";
            // std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            // std::cout << " rhoi " << rhoi << " \n";
            // int xxx;
            // std::cin >> xxx;
        }

        // if (parameters->getSplitDomain() == 0) {
        //     if (i < n_iterations/3) {
        //         if (solver_type == 0) {
        //             bicgstab_solver->set_eps(0.001);
        //             bicgstab_solver->set_newton_eps(0.1);
        //             solver->set_eps(0.001);
        //             solver->set_newton_eps(0.1);                    
        //         }
        //     }
        //     else {
        //         // continue with default settings
        //     }
        // }

        solver->solve(*epot, scharge_ave);
        efield->recalculate();

        if (parameters->getOutputVTKEnabled() && parameters->getOutputVTKExportSimulationState()) {
            const std::string& vtk_folder = fileManager->getVTKFolder();
            std::string iter_tag = fileManager->getFileTag() + "_it" + std::to_string(i+1);

            // Save epot and scharge_ave as VTK
            geometryManager->exportPotentialToVTK(*epot, vtk_folder + iter_tag);
            geometryManager->exportSpacechargeToVTK(scharge_ave, vtk_folder + iter_tag);
        }

        if (i > 0) {
            cout << " Reset trajectories...\n";
            pdb->reset_trajectories();
            cout << " Done!\n";
        }
        
        bool debugprint = false;
        if (debug && i == n_iterations-1) {
            debugprint = true;
        }	
        pdb->set_relativistic(true);
        // pdb->set_accuracy(1e-9, 1e-9);
        pdb->set_max_steps(10000);
        
        const vector<double> periodicity = parameters->getPeriodicityBounds();
        vector<double> periodicity_empty; // no periodicity considered for stripping callbacks
        
        // Initialize stripping callbacks if stripping is enabled
        THCallback_strip* thcstr = nullptr;
        THCallback_secondaries* thcsec = nullptr;

        if (parameters->getIncludeStripping() > 0) {
            // Create stripping callback objects if not already created
            if (thcstr == nullptr) {
                static double mass = parameters->getMIons();
                string density_profile = parameters->getStrippingDensityProfile();
                thcstr = new THCallback_strip(debugprint, pdb, mass, periodicity_empty, density_profile);
            }
            if (thcsec == nullptr) {
                static double mass = parameters->getMIons();
                string density_profile = parameters->getStrippingDensityProfile();
                double secondary_z_min = parameters->getStrippingMinimumZ();
                thcsec = new THCallback_secondaries(pdb, mass, periodicity, density_profile,
                                                    secondary_z_min);
            }
            
            // Set appropriate callback based on iteration and stripping mode
            if (i < n_iterations - 1) {
                pdb->set_trajectory_handler_callback(thcstr);
                if (debug) logfile << "DEBUG: Using stripping callback (iteration " << i+1 << ")" << endl;
            } else if (parameters->getIncludeStripping() > 1) {
                pdb->set_trajectory_handler_callback(thcsec);
                if (debug) logfile << "DEBUG: Using secondaries callback (final iteration " << i+1 << ")" << endl;
            } else {
                pdb->set_trajectory_handler_callback(thcstr);
                if (debug) logfile << "DEBUG: Using stripping callback (final iteration " << i+1 << ")" << endl;
            }
        }
        
        // Initialize surface collision callback (EAMCC) - only on last iteration
        static THCallback_surf_EAMCC* thc_surf_static = nullptr;
        if (parameters->getIncludeSurfaceCollisions() > 0) {
            if (i == n_iterations - 1) {
                // Last iteration: enable surface collisions
                if (thc_surf_static == nullptr) {
                    double mass = parameters->getMIons();
                    bool debug_surf = debugprint || parameters->getSurfaceCollisionsDebug();
                    Geometry* geom = geometryManager->getGeometry();
                    if (geom) {
                        thc_surf_static = new THCallback_surf_EAMCC(
                            *geom,
                            pdb,
                            mass,
                            debug_surf,
                            parameters->getSurfaceCollisionsMinimumZ());
                        pdb->set_trajectory_end_callback(thc_surf_static);
                        if (debug) logfile << "DEBUG: Surface collision callback (EAMCC) enabled on last iteration" << endl;
                    }
                } else {
                    // Reuse existing callback
                    pdb->set_trajectory_end_callback(thc_surf_static);
                    if (debug) logfile << "DEBUG: Surface collision callback (EAMCC) enabled on last iteration" << endl;
                }
            } else {
                // Not last iteration: disable surface collisions
                pdb->set_trajectory_end_callback(nullptr);
                if (debug) logfile << "DEBUG: Surface collision callback disabled (iteration " << i+1 << ")" << endl;
            }
        } else {
            // Surface collisions disabled in parameters
            pdb->set_trajectory_end_callback(nullptr);
        }
    
        
        try {
            ibsimu.message(1) << "Starting particle trajectory calculation..." << endl;
            
            // Check particle database size before iteration
            size_t initial_particles = pdb->size();
            ibsimu.message(1) << "Initial particle count: " << initial_particles << endl;
            
            // Monitor memory usage during secondary particle generation
            if (parameters->getIncludeStripping() > 1) {
                ibsimu.message(1) << "Secondary particle generation enabled - monitoring for instabilities..." << endl;
            }
            
            pdb->iterate_trajectories(*scharge, *efield, *bfield);
            
            // Check particle database size after iteration
            size_t final_particles = pdb->size();
            ibsimu.message(1) << "Final particle count: " << final_particles << endl;
            if (final_particles > initial_particles) {
                ibsimu.message(1) << "Generated " << (final_particles - initial_particles) << " secondary particles" << endl;
            }
            
            ibsimu.message(1) << "Particle trajectory calculation completed successfully." << endl;
        } catch (Error& e) {
            ibsimu.message(1) << "ERROR: Particle trajectory calculation failed" << endl;
            ibsimu.message(1) << "Error message: " << e.get_error_message() << endl;
            if (parameters->getIncludeStripping() > 1) {
                ibsimu.message(1) << "This crash occurred with secondary particle generation enabled." << endl;
                ibsimu.message(1) << "Possible causes:" << endl;
                ibsimu.message(1) << "  1. Memory overflow from too many secondary particles" << endl;
                ibsimu.message(1) << "  2. Invalid particle parameters in secondary generation" << endl;
                ibsimu.message(1) << "  3. Collision probability calculation issues" << endl;
                ibsimu.message(1) << "  4. Cross-section function returning NaN or infinity" << endl;
                ibsimu.message(1) << "Suggested fixes:" << endl;
                ibsimu.message(1) << "  1. Disable secondaries: set INCLUDE_STRIPPING = 1" << endl;
                ibsimu.message(1) << "  2. Reduce particle density or collision probabilities" << endl;
                ibsimu.message(1) << "  3. Check energy threshold conditions" << endl;
            } else {
                ibsimu.message(1) << "This usually indicates numerical instabilities in the second iteration." << endl;
                ibsimu.message(1) << "Possible solutions:" << endl;
                ibsimu.message(1) << "  1. Reduce mesh size for better stability" << endl;
                ibsimu.message(1) << "  2. Adjust space-charge averaging coefficient (ALPHA)" << endl;
                ibsimu.message(1) << "  3. Modify solver tolerances" << endl;
                ibsimu.message(1) << "  4. Use fewer iterations initially" << endl;
            }
            throw; // Re-throw to allow proper cleanup
        } catch (const std::exception& e) {
            ibsimu.message(1) << "ERROR: Standard exception in particle trajectory calculation: " << e.what() << endl;
            if (parameters->getIncludeStripping() > 1) {
                ibsimu.message(1) << "Secondary particle generation may have caused memory issues" << endl;
            }
            throw;
        } catch (...) {
            ibsimu.message(1) << "ERROR: Unknown exception in particle trajectory calculation" << endl;
            if (parameters->getIncludeStripping() > 1) {
                ibsimu.message(1) << "Unknown crash with secondary particle generation enabled" << endl;
            }
            throw;
        }
        
        rho_tot = pdb->get_rhosum();

        // Space charge averaging: scharge_ave = (scharge_old + alpha*scharge_new) / (1+alpha)
        // scharge is NOT modified in-place so the saved file contains the true last-iteration scharge.
        if (i == 0) {
            scharge_ave = *scharge;
        } else {
            double coef = parameters->getAlphaCoeff();
            MeshScalarField scharge_scaled(*scharge);
            scharge_scaled *= coef;
            scharge_ave += scharge_scaled;
            scharge_ave *= (1.0/(1.0+coef));
        }
        
        ibsimu.message(1) << " DONE!" << endl;

        bool saveITperf = parameters->getOutputSummaryEnabled();
        if (saveITperf) {
            // Update field manager with current fields
            fieldManager->setPotential(epot);
            fieldManager->setSpacecharge(scharge);
            fieldManager->setElectric(efield);
            
            vector<double> zlocitsumEG, zlocitsumOUT;
            Vec3D lastpt = geometry->max();
            zlocitsumEG.push_back(0.009);
            zlocitsumOUT.push_back(lastpt[2] - parameters->getMeshSize());

            if (debug) logfile << "DEBUG:  Iteration summary at z=" << zlocitsumEG[0] << "\n" << flush;
            if (debug) logfile << "DEBUG:  Iteration summary at z=" << zlocitsumOUT[0] << "\n" << flush;

            string iterationsfileEG = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_it_EG.txt";
            string iterationsfileOUT = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_it_OUT.txt";

            if (i==0) {
                diagnostic_data_alongZ(zlocitsumEG, i, iterationsfileEG, false);
                diagnostic_data_alongZ(zlocitsumOUT, i, iterationsfileOUT, false);
            }
            else { 
                diagnostic_data_alongZ(zlocitsumEG, i, iterationsfileEG, true);
                diagnostic_data_alongZ(zlocitsumOUT, i, iterationsfileOUT, true);
            }
            if (debug) logfile << "DEBUG:  Iteration summaries saved!\n" << flush;
        }

        if (debug) logfile << "DEBUG:  Iteration " << i+1 << "/" << n_iterations << " completed. Continuing...\n" << flush;
    }

    // Update managers with final fields
    fieldManager->setPotential(epot);
    fieldManager->setSpacecharge(scharge);
    fieldManager->setElectric(efield);

    if (parameters->getOutputDataEnabled()) {
        const string epotfile = fileManager->getDataFolder() + fileManager->getFileTag() + "_epot.dat";
        const string pdbfile = fileManager->getDataFolder() + fileManager->getFileTag() + "_pdb.dat";
        const string schargefile = fileManager->getDataFolder() + fileManager->getFileTag() + "_scharge.dat";
        const string bfieldfile = fileManager->getDataFolder() + fileManager->getFileTag() + "_bfield.dat";

        epot->save(epotfile);
        pdb->save(pdbfile);
        scharge->save(schargefile);

        // Only save magnetic field if it exists and is not the zero field we created.
        if (bfield != nullptr && zero_bfield == nullptr) {
            bfield->save(bfieldfile);
        } else {
            ofstream bfield_empty(bfieldfile);
            bfield_empty.close();
        }

        ibsimu.message(1) << " Electrostatic potential saved" << endl;
        ibsimu.message(1) << " Space-charge distribution saved" << endl;
        ibsimu.message(1) << " Electric field saved" << endl;
        ibsimu.message(1) << " Magnetic field saved" << endl;
    } else {
        ibsimu.message(1) << "Data output disabled by outputs.data settings; skipping .dat file writes" << endl;
    }

    // Cleanup zero magnetic field if we created it
    if (zero_bfield != nullptr) {
        delete zero_bfield;
        zero_bfield = nullptr;
    }

    // Cleanup dynamically allocated solver
    if (solver != nullptr) {
        delete solver;
        solver = nullptr;
    }
    ibsimu.message(1) << "Simulation completed!" << endl;
}

bool ManageSimulation::check_EGext(double oriJ, double& extsimJ) {
    return particleManager->checkEGExtractedCurrent(oriJ, extsimJ, *parameters);
}

bool ManageSimulation::load_simulation() {
    if (!parameters->getOutputDataEnabled()) {
        ibsimu.message(1) << "Loading from .dat files is disabled because outputs.data.enabled is false" << endl;
        return false;
    }

    ibsimu.message(1) << "Loading simulation from .dat files..." << endl;
    
    // Get file paths
    string data_folder = fileManager->getDataFolder();
    string file_tag = fileManager->getFileTag();
    string epotfile = data_folder + file_tag + "_1_epot.dat";
    string pdbfile = data_folder + file_tag + "_1_pdb.dat";
    string schargefile = data_folder + file_tag + "_1_scharge.dat";
    string bfieldfile = data_folder + file_tag + "_1_bfield.dat";
    
    ibsimu.message(1) << "Data folder: " << data_folder << endl;
    ibsimu.message(1) << "File tag: " << file_tag << endl;
    
    // Check if essential files exist
    ifstream test_epot(epotfile);
    ifstream test_pdb(pdbfile);
    
    if (!test_epot.good()) {
        ibsimu.message(1) << "Error: Potential file not found: " << epotfile << endl;
        test_epot.close();
        return false;
    }
    if (!test_pdb.good()) {
        ibsimu.message(1) << "Error: Particle database file not found: " << pdbfile << endl;
        test_pdb.close();
        return false;
    }
    test_epot.close();
    test_pdb.close();
    
    try {
        // First load geometry (should already be created)
        Geometry* geometry = geometryManager->getGeometry();
        if (!geometry) {
            ibsimu.message(1) << "Error: No geometry available. Create geometry first." << endl;
            return false;
        }
        ibsimu.message(1) << "Using existing geometry" << endl;
        
        // Load electrostatic potential using constructor with geometry
        ibsimu.message(1) << "Loading electrostatic potential from: " << epotfile << endl;
        std::ifstream is_epot(epotfile);
        if (!is_epot.good()) {
            throw Error(ERROR_LOCATION, "couldn't open file '" + epotfile + "'");
        }
        EpotField* epot = new EpotField(is_epot, *geometry);
        is_epot.close();
        fieldManager->setPotential(epot);
        ibsimu.message(1) << "Electrostatic potential loaded successfully" << endl;
        
        // Load space charge if available
        ifstream test_scharge(schargefile);
        if (test_scharge.good()) {
            test_scharge.close();
            ibsimu.message(1) << "Loading space charge from: " << schargefile << endl;
            std::ifstream is_scharge(schargefile);
            if (!is_scharge.good()) {
                throw Error(ERROR_LOCATION, "couldn't open file '" + schargefile + "'");
            }
            MeshScalarField* scharge = new MeshScalarField(is_scharge);
            is_scharge.close();
            fieldManager->setSpacecharge(scharge);
            ibsimu.message(1) << "Space charge loaded successfully" << endl;
        } else {
            ibsimu.message(1) << "Space charge file not found, creating zero field" << endl;
            MeshScalarField* scharge = new MeshScalarField(*geometry);
            fieldManager->setSpacecharge(scharge);
        }
        
        // Load magnetic field if available
        ifstream test_bfield(bfieldfile);
        if (test_bfield.good()) {
            test_bfield.close();
            ibsimu.message(1) << "Loading magnetic field from: " << bfieldfile << endl;
            std::ifstream is_bfield(bfieldfile);
            if (!is_bfield.good()) {
                throw Error(ERROR_LOCATION, "couldn't open file '" + bfieldfile + "'");
            }
            MeshVectorField* bfield = new MeshVectorField(is_bfield);
            is_bfield.close();
            fieldManager->setMagnetic(bfield);
            ibsimu.message(1) << "Magnetic field loaded successfully" << endl;
        } else {
            ibsimu.message(1) << "Magnetic field file not found, creating zero field" << endl;
            bool fout[3] = {true, true, true};
            MeshVectorField* bfield = new MeshVectorField(*geometry, fout);
            fieldManager->setMagnetic(bfield);
        }
        
        // Create electric field from potential
        ibsimu.message(1) << "Creating electric field from potential..." << endl;
        EpotEfield* efield = new EpotEfield(*epot);
        field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                         FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                         FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
        efield->set_extrapolation(efldextrpl);
        fieldManager->setElectric(efield);
        ibsimu.message(1) << "Electric field created successfully" << endl;
        
        // Load particle database using constructor with geometry
        ibsimu.message(1) << "Loading particle database from: " << pdbfile << endl;
        std::ifstream is_pdb(pdbfile);
        if (!is_pdb.good()) {
            throw Error(ERROR_LOCATION, "couldn't open file '" + pdbfile + "'");
        }
        ParticleDataBase3D* pdb = new ParticleDataBase3D(is_pdb, *geometry);
        is_pdb.close();
        particleManager->setParticles(pdb);
        ibsimu.message(1) << "Particle database loaded successfully with " << pdb->size() << " particles" << endl;
        
        ibsimu.message(1) << "All simulation data loaded successfully!" << endl;
        return true;
        
    } catch (const Error& e) {
        ibsimu.message(1) << "IBSimu Error during loading" << endl;
        return false;
    } catch (const exception& e) {
        ibsimu.message(1) << "Standard exception during loading: " << e.what() << endl;
        return false;
    } catch (...) {
        ibsimu.message(1) << "Unknown exception during loading" << endl;
        return false;
    }
}

bool ManageSimulation::trace_particles_with_loaded_fields(bool use_stripping) {
    ibsimu.message(1) << "Starting particle tracing with loaded fields..." << endl;
    
    try {
        // Verify all required fields are loaded
        Geometry* geometry = geometryManager->getGeometry();
        EpotField* epot = fieldManager->getPotential();
        EpotEfield* efield = fieldManager->getElectric();
        MeshVectorField* bfield = fieldManager->getMagnetic();
        MeshScalarField* scharge = fieldManager->getSpacecharge();
        ParticleDataBase3D* pdb = particleManager->getParticles();
        
        if (!geometry || !epot || !efield || !bfield || !scharge || !pdb) {
            ibsimu.message(1) << "Error: Missing required fields or geometry for particle tracing" << endl;
            ibsimu.message(1) << "Geometry: " << (geometry ? "OK" : "MISSING") << endl;
            ibsimu.message(1) << "Potential: " << (epot ? "OK" : "MISSING") << endl;
            ibsimu.message(1) << "Electric: " << (efield ? "OK" : "MISSING") << endl;
            ibsimu.message(1) << "Magnetic: " << (bfield ? "OK" : "MISSING") << endl;
            ibsimu.message(1) << "Space charge: " << (scharge ? "OK" : "MISSING") << endl;
            ibsimu.message(1) << "Particles: " << (pdb ? "OK" : "MISSING") << endl;
            return false;
        }
        
        ibsimu.message(1) << "All fields verified. Starting particle tracing..." << endl;
        ibsimu.message(1) << "Number of particles to trace: " << pdb->size() << endl;
        
        // Set up magnetic field suppression
        double usup = 10.0; 
        NPlasmaBfieldSuppression psup(*epot, usup);
        pdb->set_bfield_suppression(&psup);
        
        // Set relativistic mode
        pdb->set_relativistic(true);
        
        // Initialize stripping callbacks if requested
        THCallback_strip* thcstr = nullptr;
        THCallback_secondaries* thcsec = nullptr;
        
        if (use_stripping && parameters->getIncludeStripping() > 0) {
            ibsimu.message(1) << "Setting up stripping callbacks..." << endl;
            
            double mass = parameters->getMIons();
            string density_profile = parameters->getStrippingDensityProfile();
            vector<double> periodicity_empty; // no periodicity for main stripping
            
            thcstr = new THCallback_strip(false, pdb, mass, periodicity_empty, density_profile);
            
            // if (parameters->getIncludeStripping() > 1) {
            if (true) { // Always create secondaries callback for tracing if stripping is enabled
                const vector<double> periodicity = parameters->getPeriodicityBounds();
                double secondary_z_min = parameters->getStrippingMinimumZ();
                thcsec = new THCallback_secondaries(pdb, mass, periodicity, density_profile,
                                                    secondary_z_min);
                pdb->set_trajectory_handler_callback(thcsec);
                ibsimu.message(1) << "Using secondaries stripping callback" << endl;
                
                // Add warning about potential memory issues
                ibsimu.message(1) << "WARNING: Secondary particle generation enabled" << endl;
                ibsimu.message(1) << "This may cause memory issues or crashes with large particle numbers" << endl;
                ibsimu.message(1) << "Current particle count: " << pdb->size() << endl;
                
            } else {
                pdb->set_trajectory_handler_callback(thcstr);
                ibsimu.message(1) << "Using basic stripping callback" << endl;
            }
        } else {
            ibsimu.message(1) << "Stripping disabled for particle tracing" << endl;
        }
        
        // Initialize surface collision callback (EAMCC secondary particle generation)
        THCallback_surf_EAMCC* thc_surf = nullptr;
        
        if (parameters->getIncludeSurfaceCollisions() > 0) {
            ibsimu.message(1) << "Setting up surface collision callback (EAMCC)..." << endl;
            
            double mass = parameters->getMIons();
            bool debug_surf = debug || parameters->getSurfaceCollisionsDebug();
            
            thc_surf = new THCallback_surf_EAMCC(
                *geometry,
                pdb,
                mass,
                debug_surf,
                parameters->getSurfaceCollisionsMinimumZ());
            pdb->set_trajectory_end_callback(thc_surf);
            
            ibsimu.message(1) << "Surface collision callback enabled (EAMCC model)" << endl;
            ibsimu.message(1) << "  - Electron backscattering and secondary emission" << endl;
            ibsimu.message(1) << "  - Ion backscattering and secondary electron emission" << endl;
            ibsimu.message(1) << "  - Surface-generated particles use generation offset +101" << endl;
        } else {
            ibsimu.message(1) << "Surface collisions disabled" << endl;
        }
        
        // Reset trajectories before tracing
        pdb->reset_trajectories();
        ibsimu.message(1) << "Particle trajectories reset" << endl;
        
        // Perform particle tracing
        ibsimu.message(1) << "Iterating particle trajectories..." << endl;
        pdb->iterate_trajectories(*scharge, *efield, *bfield);
        ibsimu.message(1) << "Particle trajectory iteration completed successfully!" << endl;
        
        if (parameters->getOutputDataEnabled()) {
            const string pdb_traced_file = fileManager->getDataFolder() + fileManager->getFileTag() + "_traced_pdb.dat";
            pdb->save(pdb_traced_file);
            ibsimu.message(1) << "Traced particle database saved to: " << pdb_traced_file << endl;
        }
        
        // Export traced particle trajectories
        if (parameters->getOutputVTKEnabled() && parameters->getOutputVTKExportTracedParticles() &&
            particleManager && particleManager->getParticles() && particleManager->getParticles()->size() > 0) {
            const string vtk_base = fileManager->getVTKFolder() + fileManager->getFileTag() + "_traced";
            particleManager->exportTrajectoriesToVTK(vtk_base);
            ibsimu.message(1) << "Traced trajectories exported to VTK: " << vtk_base << "_trajectories.vtk" << endl;
        }
        
        // Cleanup stripping callbacks
        if (thcstr) delete thcstr;
        if (thcsec) delete thcsec;
        
        ibsimu.message(1) << "Particle tracing with loaded fields completed successfully!" << endl;
        return true;
        
    } catch (const Error& e) {
        ibsimu.message(1) << "IBSimu Error during particle tracing" << endl;
        return false;
    } catch (const exception& e) {
        ibsimu.message(1) << "Standard exception during particle tracing: " << e.what() << endl;
        return false;
    } catch (...) {
        ibsimu.message(1) << "Unknown exception during particle tracing" << endl;
        return false;
    }
}

void ManageSimulation::diagnostic_data_alongZ(const vector<double>& diagzpos, int iteration, 
                                             const std::string& outfile, bool append_at_end) {
    // Use the enhanced diagnostics function with field calculations
    diagnosticsManager->generateDiagnosticData(diagzpos, iteration, outfile, append_at_end,
                                              particleManager->getParticles(),
                                              geometryManager->getGeometry(),
                                              fieldManager->getPotential(),
                                              fieldManager->getMagnetic());
}

void ManageSimulation::diagnostic_data_alongZ(const vector<double>& diagzpos, int iteration, 
                                             const std::string& outfile) {
    diagnostic_data_alongZ(diagzpos, iteration, outfile, false);
}

void ManageSimulation::plot_simulation(int argc, char **argv) {
    if (!parameters->getOutputPlotsEnabled()) {
        ibsimu.message(1) << "Plot generation disabled by outputs.plots settings" << endl;
        return;
    }

    diagnosticsManager->createPlots(argc, argv, 
                                   geometryManager->getGeometry(),
                                   fieldManager->getPotential(),
                                   fieldManager->getMagnetic(),
                                   fieldManager->getElectric(),
                                   fieldManager->getSpacecharge(),
                                   particleManager->getParticles(),
                                   fileManager->getPlotFolder(),
                                   fileManager->getFileTag());
}

void ManageSimulation::plot_simulation(int argc, char **argv, particle_kind pk) {
    if (!parameters->getOutputPlotsEnabled()) {
        ibsimu.message(1) << "Plot generation disabled by outputs.plots settings" << endl;
        return;
    }

    diagnosticsManager->createPlots(argc, argv,
                                   geometryManager->getGeometry(),
                                   fieldManager->getPotential(),
                                   fieldManager->getMagnetic(),
                                   fieldManager->getElectric(),
                                   fieldManager->getSpacecharge(),
                                   particleManager->getParticles(),
                                   particleManager->getParticleSpecies(),
                                   pk,
                                   fileManager->getPlotFolder(),
                                   fileManager->getFileTag());
}

void ManageSimulation::analysis(double zlocsummary) {
    if (!parameters->getOutputSummaryEnabled()) {
        ibsimu.message(1) << "Summary output disabled by outputs.summary settings" << endl;
        return;
    }

    // Use enhanced analysis with field calculations
    diagnosticsManager->performAnalysis(particleManager->getParticles(), *parameters,
                                       geometryManager->getGeometry(),
                                       fieldManager->getPotential(),
                                       fieldManager->getMagnetic(),zlocsummary,
                                       particleManager->getParticleSpecies(),
                                       true, fileManager.get(),
                                       geometryManager->getZGrids());
}

void ManageSimulation::fill_particle_dbs() {
    particleManager->fillParticleDatabases(*parameters);
}

void ManageSimulation::print_traj_to_txt(const std::string& filename) {
    diagnosticsManager->printTrajectoryData(filename, particleManager->getParticles());
}

std::string ManageSimulation::save_emitter(const std::string& emitname, double zloc) {
    if (!parameters->getOutputSummaryEnabled()) {
        ibsimu.message(1) << "Emitter export disabled by outputs.summary settings" << endl;
        return std::string();
    }

    return particleManager->saveEmitter(fileManager->getFileTag()+emitname, zloc, fileManager->getOutputSummaryFolder());
}

void ManageSimulation::particles_end_location() {
    particleManager->analyzeParticleEndLocations(*parameters);
}

std::vector<double> ManageSimulation::analyze_grid_power_loads(particle_kind pk) {
    if (!particleManager->getParticles() || !geometryManager->getGeometry()) {
        ibsimu.message(1) << "Error: Missing particle database or geometry for grid power analysis" << endl;
        return std::vector<double>();
    }
    
    // Get grid positions and mesh size from parameters
    std::vector<double> zgrids = geometryManager->getZGrids();
    double mesh_size = parameters->getMeshSize();
    double ionmass = parameters->getMIons();
    
    // Use default grid positions if none available
    if (zgrids.empty()) {
        zgrids = {0.0, 0.009, 0.015, 0.032, 0.120, 0.137, 0.225, 0.242, 0.330, 0.347, 0.435, 0.452, 0.540, 0.557};
        ibsimu.message(1) << "Using default MTF grid positions for power analysis" << endl;
    }
    
    ibsimu.message(1) << "Analyzing grid power loads..." << endl;
    
    // Call the diagnostics manager to perform the analysis
    std::vector<double> power_results = diagnosticsManager->analyzeGridPowerLoads(
        particleManager->getParticles(),
        zgrids,
        mesh_size,
        geometryManager->getGeometry(),
        ionmass,
        fileManager->getOutputSummaryFolder(),
        fileManager->getFileTag(),
        pk
    );
    
    ibsimu.message(1) << "Grid power analysis completed. Found " << power_results.size() << " power values." << endl;
    
    return power_results;
}

void ManageSimulation::create_simulation_summary(double zlocsummary, const std::string& summary_file_tag, 
                                                bool append_at_end) {
    // Calculate extracted current density at EG exit using ParticleManager
    double oriJ = parameters->getJIon();  // Original input current density
    double extsimJ = 0.0;  // Will be filled by checkEGExtractedCurrent
    
    // Get the extracted current density
    bool eg_check = particleManager->checkEGExtractedCurrent(oriJ, extsimJ, *parameters);
    
    if (debug) {
        logfile << "DEBUG: EG extracted current check: " << (eg_check ? "PASS" : "FAIL") << endl;
        logfile << "DEBUG: Target J=" << oriJ << " A/m², Extracted J=" << extsimJ << " A/m²" << endl;
    }
    
    // Use the enhanced version with field calculations and extracted current density
    diagnosticsManager->createSimulationSummary(zlocsummary, summary_file_tag, append_at_end,
                                               particleManager->getParticles(),
                                               fileManager->getOutputSummaryFolder(),
                                               geometryManager->getGeometry(),
                                               fieldManager->getPotential(),
                                               fieldManager->getMagnetic(),
                                               extsimJ);  // Pass extracted current density
}

void ManageSimulation::create_simulation_summary(double zlocsummary, const std::string& summary_file_tag, 
                                                bool append_at_end, particle_kind pk) {
    diagnosticsManager->createSimulationSummary(zlocsummary, summary_file_tag, append_at_end,
                                               particleManager->getParticleSpecies(), pk,
                                               fileManager->getOutputSummaryFolder());
}

// Individual simulation summary for scan
void ManageSimulation::create_individual_simulation_summary(int scan_index, const std::string& simulation_tag,
                                         double zlocsummary, particle_kind pk) {
    // Calculate extracted current density at EG exit using ParticleManager
    double oriJ = parameters->getJIon();  // Original input current density
    double extsimJ = 0.0;  // Will be filled by checkEGExtractedCurrent
    
    // Get the extracted current density
    bool eg_check = particleManager->checkEGExtractedCurrent(oriJ, extsimJ, *parameters);
    
    if (debug) {
        logfile << "DEBUG: Individual summary EG extracted current check: " << (eg_check ? "PASS" : "FAIL") << endl;
        logfile << "DEBUG: Target J=" << oriJ << " A/m², Extracted J=" << extsimJ << " A/m²" << endl;
    }
    
    diagnosticsManager->createIndividualSimulationSummary(scan_index, zlocsummary,
                                                         particleManager->getParticles(),
                                                         extsimJ,
                                                         fileManager.get(),
                                                         pk);
}

// Add to scan-level beam properties summary
void ManageSimulation::add_to_scan_beam_properties_summary(int scan_index, const std::string& simulation_tag,
                                        const std::string& scan_folder, const std::string& scan_file_tag,
                                        double zlocsummary, particle_kind pk) {
    // Calculate extracted current density at EG exit using ParticleManager
    double oriJ = parameters->getJIon();  // Original input current density
    double extsimJ = 0.0;  // Will be filled by checkEGExtractedCurrent
    const ParticleDataBase3D* beam_property_particles = particleManager->getParticles();
    const std::vector<ParticleDataBase3D*>& particle_species = particleManager->getParticleSpecies();
    int negion_index = get_particle_int(PARTICLE_HM);
    if (negion_index >= 0 && negion_index < static_cast<int>(particle_species.size()) && particle_species[negion_index]) {
        beam_property_particles = particle_species[negion_index];
    }
    
    // Get the extracted current density
    particleManager->checkEGExtractedCurrent(oriJ, extsimJ, *parameters);
    
    diagnosticsManager->addToScanBeamPropertiesSummary(scan_index, simulation_tag, scan_folder, scan_file_tag,
                                                      zlocsummary, beam_property_particles,
                                                      extsimJ, fieldManager->getMagnetic(),
                                                      fieldManager->getPotential(), pk);
}

// Add to scan-level grid power summary
void ManageSimulation::add_to_scan_grid_power_summary(int scan_index, const std::string& simulation_tag,
                                     const std::string& scan_folder, const std::string& scan_file_tag,
                                     particle_kind pk) {
    diagnosticsManager->addToScanGridPowerSummary(scan_index, simulation_tag, scan_folder, scan_file_tag,
                                                 particleManager->getParticles(), 
                                                 geometryManager->getGeometry(), pk);
}

void ManageSimulation::set_file_tag(const std::string& filetag) {
    fileManager->setFileTag(filetag);
}

// PowerStruct implementation (kept for compatibility)
void PowerStruct::add(const Particle3D &pp) {
    Vec3D x = pp.location();
    Vec3D vel = pp.velocity();
    size_t ps = identify_particle_species(pp.m(), pp.q(), ionmass);
    bool addpart = true;
    if (ps == PARTICLE_HM && vel[2] < 0) addpart = false;

    if (addpart) {
        xdata.push_back(x[0]);
        ydata.push_back(x[1]);
        zdata.push_back(x[2]);

        kinddata.push_back(static_cast<int>(ps));
        vxdata.push_back(vel[0]);
        vydata.push_back(vel[1]);
        vzdata.push_back(vel[2]);
        curdata.push_back(pp.IQ());
        mdata.push_back(pp.m());
        double charge = pp.q();
        if (charge == 0) charge = CHARGE_E;
        double V = 0.5 * pp.m() * vel.ssqr() / charge;
        wdata.push_back(pp.IQ() * V);
        qdata.push_back(pp.q());
        gendata.push_back(pp.gen());
    }
}

void PowerStruct::calculate_total_power() {
    total_power = 0.;
    total_current = 0.;

    total_power_perspecies.assign(6, 0.0);
    total_current_perspecies.assign(6, 0.0);

    total_power_pergen.assign(6, 0.0);
    total_current_pergen.assign(6, 0.0);

    for (size_t ii = 0; ii < wdata.size(); ii++) {
        total_power += wdata[ii];
        total_current += curdata[ii];
        if (kinddata[ii] >= 0 && kinddata[ii] < static_cast<int>(total_power_perspecies.size())) {
            total_power_perspecies[kinddata[ii]] += wdata[ii];
            total_current_perspecies[kinddata[ii]] += curdata[ii];
        }
        if (abs(gendata[ii] % 100) < 5) {
            total_power_pergen[gendata[ii]] += wdata[ii];
            total_current_pergen[gendata[ii]] += curdata[ii];
        } else {
            total_power_pergen[5] += wdata[ii];
            total_current_pergen[5] += curdata[ii];
        }
    }
}

void PowerStruct::print(const std::string &filename) {
    ofstream ostr(filename.c_str());
    if (!ostr.is_open()) {
        cerr << "Error: Could not open file " << filename << " for writing\n";
        return;
    }
    
    ostr << "x\ty\tz\tvx\tvy\tvz\tm\tpower\tcurr\tcharge\tgen\tkind\n";
    for (size_t i = 0; i < xdata.size(); i++) {
        ostr << setw(12) << xdata[i] << "\t"
             << setw(12) << ydata[i] << "\t"
             << setw(12) << zdata[i] << "\t"
             << setw(12) << vxdata[i] << "\t"
             << setw(12) << vydata[i] << "\t"
             << setw(12) << vzdata[i] << "\t"
             << setw(12) << mdata[i] << "\t"
             << setw(12) << wdata[i] << "\t"
             << setw(12) << curdata[i] << "\t"
             << setw(12) << qdata[i] << "\t"
             << setw(12) << gendata[i] << "\t"
             << setw(12) << kinddata[i] << "\n";
    }
    ostr << "\n";
    ostr.close();
}

void ManageSimulation::export_for_paraview(const std::string& base_filename) {
    ibsimu.message(1) << "Exporting simulation data for ParaView visualization..." << endl;
    
    // Export geometry
    if (geometryManager && geometryManager->getGeometry()) {
        // Only export solid boundaries during simulation run
        // (Other VTK export functions remain available for debugging)
        geometryManager->exportSolidsToVTK(base_filename);
    } else {
        ibsimu.message(1) << "Warning: No geometry available for export" << endl;
    }
    
    // Export potential field if available
    if (fieldManager && fieldManager->getPotential()) {
        geometryManager->exportPotentialToVTK(*fieldManager->getPotential(), base_filename);
    } else {
        ibsimu.message(1) << "Warning: No potential field available for export" << endl;
    }
    
    ibsimu.message(1) << "ParaView export completed with base filename: " << base_filename << endl;
    ibsimu.message(1) << "Files created:" << endl;
    ibsimu.message(1) << "  - " << base_filename << "_solids.vtk (solid boundaries)" << endl;
    ibsimu.message(1) << "  - " << base_filename << "_potential.vtk (potential field)" << endl;
}

void ManageSimulation::export_geometry_to_vtk(const std::string& base_filename) {
    if (!parameters->getOutputVTKEnabled() || !parameters->getOutputVTKExportGeometry()) {
        ibsimu.message(1) << "Geometry VTK export disabled by outputs.vtk settings" << endl;
        return;
    }

    // Create VTK directory if it doesn't exist
    string vtk_dir = base_filename.substr(0, base_filename.find_last_of('/'));
    string mkdir_cmd = "mkdir -p " + vtk_dir;
    system(mkdir_cmd.c_str());
    
    ibsimu.message(1) << "Exporting geometry to VTK format at simulation start..." << endl;
    ibsimu.message(1) << "VTK directory: " << vtk_dir << endl;
    
    // Export geometry if available
    if (geometryManager && geometryManager->getGeometry()) {
        // Only export solid boundaries during simulation run
        // (Other VTK export functions remain available for debugging)
        geometryManager->exportSolidsToVTK(base_filename);
        
        ibsimu.message(1) << "Geometry VTK export completed!" << endl;
        ibsimu.message(1) << "Files created in " << vtk_dir << ":" << endl;
        ibsimu.message(1) << "  - " << base_filename.substr(base_filename.find_last_of('/')+1) << "_solids.vtk (solid boundaries)" << endl;
    } else {
        ibsimu.message(1) << "Warning: No geometry available for VTK export" << endl;
    }
}

void ManageSimulation::save_results_to_vtk() {
    if (!parameters->getOutputVTKEnabled()) {
        ibsimu.message(1) << "VTK export disabled by outputs.vtk settings" << endl;
        return;
    }

    if (!parameters->getOutputVTKExportSimulationState() &&
        !parameters->getOutputVTKExportTracedParticles() &&
        !parameters->getOutputVTKExportGeometry()) {
        ibsimu.message(1) << "All end-of-run VTK exports disabled by outputs.vtk settings" << endl;
        return;
    }

    // Export complete simulation results to VTK format for ParaView visualization
    string vtk_base = fileManager->getVTKFolder() + fileManager->getFileTag() + "_simulation";

    export_simulation_results_to_vtk(vtk_base, fieldManager->getSpacecharge());

}

void ManageSimulation::export_simulation_results_to_vtk(const std::string& base_filename,
                                                        const MeshScalarField* scharge) {
    if (!geometryManager || !geometryManager->getGeometry()) {
        ibsimu.message(1) << "Warning: No geometry available for VTK export" << endl;
        return;
    }
    
    // Create VTK directory if it doesn't exist
    string vtk_dir = base_filename.substr(0, base_filename.find_last_of('/'));
    string mkdir_cmd = "mkdir -p " + vtk_dir;
    system(mkdir_cmd.c_str());
    
    ibsimu.message(1) << "Exporting complete simulation results to VTK format..." << endl;
    
    // Export space-charge field
    if (parameters->getOutputVTKExportSimulationState()) {
        if (scharge) {
            geometryManager->exportSpacechargeToVTK(*scharge, base_filename);
        } else {
            ibsimu.message(1) << "Warning: No space-charge field available for VTK export" << endl;
        }
    }
    
    // Export particle trajectories
    if (parameters->getOutputVTKExportTracedParticles() &&
        particleManager && particleManager->getParticles() && particleManager->getParticles()->size() > 0) {
        particleManager->exportTrajectoriesToVTK(base_filename);
    } else if (parameters->getOutputVTKExportTracedParticles()) {
        ibsimu.message(1) << "Warning: No particle trajectories available for VTK export" << endl;
    }
    
    // Export potential field (already available but let's add it for completeness)
    if (parameters->getOutputVTKExportSimulationState()) {
        if (fieldManager && fieldManager->getPotential()) {
            geometryManager->exportPotentialToVTK(*fieldManager->getPotential(), base_filename);
        } else {
            ibsimu.message(1) << "Warning: No potential field available for VTK export" << endl;
        }
    }
    
    // Export geometry boundaries when requested.
    if (parameters->getOutputVTKExportGeometry()) {
        geometryManager->exportSolidsToVTK(base_filename);
    }
    
    ibsimu.message(1) << "VTK export completed! Files created in " << vtk_dir << ":" << endl;
    ibsimu.message(1) << "  - " << base_filename.substr(base_filename.find_last_of('/')+1) << "_scharge.vtk (space-charge field)" << endl;
    ibsimu.message(1) << "  - " << base_filename.substr(base_filename.find_last_of('/')+1) << "_trajectories.vtk (particle trajectories)" << endl;
    ibsimu.message(1) << "  - " << base_filename.substr(base_filename.find_last_of('/')+1) << "_potential.vtk (potential field)" << endl;
    ibsimu.message(1) << "  - " << base_filename.substr(base_filename.find_last_of('/')+1) << "_solids.vtk (solid boundaries)" << endl;
}
