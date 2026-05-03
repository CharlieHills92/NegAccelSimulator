#ifndef _PARAMETERS
#define _PARAMETERS

// Simulation OUTPUT 
std::string output_name="SPIDER";  // Name (of the folder containing the results);
bool full_output = false;   // if true also outputs the potential and charge map and the particle database. Output size can be huge in multi-scans. 
bool save_pot_txt = false;  // save a txt file wiht 3d distrib of the potential (only one, must chnage the string if map for all cases during a scan are needed.


//parametri 
int numero_cicli = 8;                                 //iteration number
double h = 0.0005;                                   //Mesh cell size [m]

int b_ison = 1;                                       //=1 if B is used
const std::string bfieldfn = "./BFields/SPIDER_Rev_CESM_GGFe.txt";  // name of B field map
double bscaling=1.0; // scale the B field by a certain amount in all directions
double sc_alpha = 0.8;                                //Space charge averaging factor
//beam parameters
int ion_mass=1;                    // ion mass (amu)
double ion_charge=-1;          //ion charge
int numero_ioni = 30000;                 //number of ions
int numero_elettroni = 0;                //number of electrons (not sure if it works!)
double el_current_density = 0.0;         //Densità di corrente degli elettroni [A/m2]
double h_current_density = 100;       //(obsolete/unused in scans) Densità di corrente degli ioni [A/m2]
double beam_mean_energy = 3;                        //beam_mean_energy in z direction[eV]
double beam_Tp = 0;                                 //Beam parallel temperature [eV]
double beam_Tt = 0;                                 //Beam transverse temperature [eV]
double beam_radius = 0.012;                           //Beam initial radius in the plasma [m]
//plasma parameters
double max_plasma = 0.007;                          //First guess of plasma meniscus z coordinate [m], also used to evaluate extracted current
double Tp = 0.8;             // Temp of electrons in pos ions simulation (P. Veltri 02052019   //Temp. f positive ions; first fast protons, then thermal ions
double Rf = 1;                                      //Ratio of fast positive ions to total negative charge density
double Up = 1;                                      //Plasma potential
double usup = 10.0;                                   //Threshold value for magnetic suppression inside plasma


// Space Charge Compensation
double SCC_degree=1;
double SCC_z=0.08;




//Plotting Parameters
int part_graph = 1;                                   //if !=0 plot graph at each iterations
int traj2plot=500;			// How many trajectories will be plotted
double equip2plot[5]={1, 5, 10, 100, 1000};			// Values of equipotential to be plotted
int Nzplane=25;			// number of points where to save divergence (PG area will be automatically refined)


// param for 2D simulation with scans


double z_diagn=0.0785; // were output (divergence, etc.. ) is evaluated 9typically 5 or 10 mm after GG)


// geom

//double z_shift=0.003; //shift all geometry by x mm

const int N_aperture_x=1;               //Number apertures in  x 
const int N_aperture_y=1;               //Number apertures in  y
const int numero_fori = N_aperture_x*N_aperture_y;    //Numero fori delle griglie  tot
const int numero_griglie = 3;                         //number of electrodes

double Lx05 = 0.01; 		// half Cell size in x (m) to set periodicity
double Ly05 = 0.011;            // half Cell size in y (m) to set periodicity

double xC_fori[numero_fori] = {0.0}; // center of each aperture in x 
double yC_fori[numero_fori] = {0.0}; // center of each aperture in y 



//mesh 
double mesh_x = 0.026;                                 //Mesh x_dimension [m]
double mesh_y = 0.028;                               //Mesh y_dimension [m]
double mesh_z = 0.08;                                //Mesh z_dimension [m]



// VOltages and Current
double V0=0.0;// grid 0 voltage in Volts (tipically PG)

// grid voltages and current density used for the scans

const int nscan=1;// number of simulations 
double current_density_scan[nscan]={330}; // extracted currents, A/m2
double V1_scan[nscan]={9.8}; // grid 1 voltage (tipically EG)
double V2_scan[nscan]={100};// grid 2 voltage (tipically AG1/GG)

/*
const int nscan=1;// number of simulations 
double current_density_scan[nscan]={150}; // extracted currents, A/m2
double V1_scan[nscan]={5000}; // grid 1 voltage (tipically EG)
double V2_scan[nscan]={31000};// grid 2 voltage (tipically AG1/GG)

*/

//potenziale griglie [V] (potenziale delle griglie in ordine da quella a contatto con il plasma)
//double V_grids[numero_griglie] = {0.0,5000.0,33000.0};     // unused                
//posizione lungo l'asse (coordinata z) delle griglie [m]
//double z_grids[numero_griglie] = {0.003,0.0098,0.0153};                      
//double z_grids[numero_griglie] = {0.00,0.0073,0.0128};  // unused in 2D
//coordinate dei centri dei fori delle griglie [m]





#endif

