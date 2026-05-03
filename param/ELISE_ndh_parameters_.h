
// Simulation OUTPUT 
string output_name="ELI_1x1_EG5_j175_Tp0_Tt0_B0_AGscan_B0";  // Name (of the folder containing the results);
bool full_output = false;   // if true also outputs the potential and charge map and the particle database. Output size can be huge in multi-scans. 



//parametri 
int numero_cicli = 20;                                 //iteration number
double h = 0.0005;                                   //Mesh cell size [m]
int part_graph = 0;                                   //if !=0 output graph at each iterations
int b_ison = 0;                                       //=1 if B is used
const std::string bfieldfn = "ELISE_1beamlet_B.txt";  // name of B field map
double sc_alpha = 0.8;                                //Space charge averaging factor
//beam parameters
int ion_mass=1;                    // ion mass (amu)
double ion_charge=-1;          //ion charge
int numero_ioni = 50000;                               //number of ions
int numero_elettroni = 0;                           //number of electrons
double el_current_density = 0.0;                   //Densità di corrente degli elettroni [A/m2]
double h_current_density = 100;                      //(obsolete/unused in scans) Densità di corrente degli ioni [A/m2]
double beam_mean_energy = 3;                        //beam_mean_energy [eV]
double beam_Tp = 0.0;                                 //Beam parallel temperature [eV]
double beam_Tt = 0.0;                                 //Beam transverse temperature [eV]
double beam_radius = 0.01;                           //Beam initial radius [m]
//plasma parameters
double max_plasma = 0.007;                          //First guess of plasma meniscus z coordinate [m], also used to evaluate extracted current
double Tp = 3;             // Temp of electrons in pos ions simulation (P. Veltri 02052019   //Temp. f positive ions; first fast protons, then thermal ions
double Rf = 0.2;                                      //Ratio of fast positive ions to total negative charge density
double Up = 1;                                      //Plasma potential
double usup = 10.0;                                   //Threshold value for magnetic suppression inside plasma

// param for 2D simulation with scans


double z_diagn=0.0625; // were output (divergence, etc.. ) is evaluated


// geom

const int N_aperture_x=1;               //Number apertures in  x 
const int N_aperture_y=1;               //Number apertures in  y
const int numero_fori = N_aperture_x*N_aperture_y;    //Numero fori delle griglie  tot
const int numero_griglie = 3;                         //number of electrodes
double Lx05 = 0.01; 		// half Cell size in x (m) to set periodicity
double Ly05 = 0.01;            // half Cell size in y (m) to set periodicity

double xC_fori[numero_fori] = {0.0}; // center of each aperture in x 
double yC_fori[numero_fori] = {0.0}; // center of each aperture in y 



//mesh 
double mesh_x = 0.025;                                 //Mesh x_dimension [m]
double mesh_y = 0.025;                               //Mesh y_dimension [m]
double mesh_z = 0.063;                                //Mesh z_dimension [m]



// VOltages and Current
double V0=0.0;// grid 0 voltage in Volts (tipically PG)

// grid voltages and current density used for the scans
const int nscan=5;// number of simulations 
double current_density_scan[nscan]={175,175,175,175,175}; // extracted currents, A/m2
double V1_scan[nscan]={5000,5000,5000,5000,5000}; // grid 1 voltage (tipically EG)
double V2_scan[nscan]={29000,31000,33000,35000,37000};// grid 2 voltage (tipically AG1/GG)


//potenziale griglie [V] (potenziale delle griglie in ordine da quella a contatto con il plasma)
//double V_grids[numero_griglie] = {0.0,5000.0,33000.0};     // unused                
//posizione lungo l'asse (coordinata z) delle griglie [m]
//double z_grids[numero_griglie] = {0.003,0.0098,0.0153};                      
//double z_grids[numero_griglie] = {0.00,0.0073,0.0128};  // unused in 2D
//coordinate dei centri dei fori delle griglie [m]

