#ifndef _PARAMETERS
#define _PARAMETERS

#include "../globals.h"
#include <string>

// Simulation OUTPUT 
std::string output_name="TEST";  // Name (of the folder containing the results);
bool full_output = false;   // if true also outputs the potential and charge map and the particle database. Output size can be huge in multi-scans. 
bool save_pot_txt = false;  // save a txt file wiht 3d distrib of the potential (only one, must chnage the string if map for all cases during a scan are needed.



//parametri 
int numero_cicli = 15;                                 //iteration number
double h = 0.001;                                   //Mesh cell size [m]

int b_ison =0;                                       //=1 if B is used
const std::string bfieldfn = "./BFields/ELISE_1beamlet_B.txt";  // name of B field map
double bscaling=1.0; // scale the B field by a certain amount in all directions
double sc_alpha = 0.8;                                //Space charge averaging factor
//beam parameters
int ion_mass=1;                    // ion mass (amu)
double ion_charge=-1;          //ion charge
int numero_ioni = 50000;                               //number of ions
int numero_elettroni = 0;                           //number of electrons
double el_current_density = 0.0;                   //Densità di corrente degli elettroni [A/m2]
double h_current_density = 100;                      //(obsolete/unused in scans) Densità di corrente degli ioni [A/m2]
double beam_mean_energy = 3;                        //beam_mean_energy [eV]
double beam_Tp = 0.000;                                 //Beam parallel temperature [eV]
double beam_Tt = 1.0;                                 //Beam transverse temperature [eV]
double beam_radius = 0.012;                           //Beam initial radius [m]
//plasma parameters
double max_plasma = 0.007;                          //First guess of plasma meniscus z coordinate [m], also used to evaluate extracted current
double Tp = 0.8;             // Temp of electrons in pos ions simulation (P. Veltri 02052019   //Temp. f positive ions; first fast protons, then thermal ions
double Rf = 1;                                      //Ratio of fast positive ions to total negative charge density
double Up = 1;                                      //Plasma potential
double usup = 10.0;                                   //Threshold value for magnetic suppression inside plasma

// Space Charge Compensation
double SCC_degree=1;
double SCC_z=0.55;

//Plotting Parameters
int part_graph = 0;                                   //if !=0 output graph at each iterations
int traj2plot=500;			// How many trajectories will be plotted
double equip2plot[5]={1, 5, 10, 100, 1000};			// Values of equipotential to be plotted
int Nzplane=80;			// number of points where to plot divergence (automatic step of 1mm for the first 30 mm)


// param for 2D simulation with scans


double z_diagn=0.553; // were output (divergence, etc.. ) is evaluated (usually 10 mm after GG)


// geom

//double z_shift=0.003; //shift all geometry by x mm

// const int N_aperture_x=1;               //Number apertures in  x 
// const int N_aperture_y=1;               //Number apertures in  y
const int numero_fori = N_aperture_x*N_aperture_y;    //Numero fori delle griglie  tot

# define _Many_electrodes  // to be defined if more than 3 electrodes are used (to model mitica/HNB)
const int numero_griglie = 7;                         //number of electrodes
double Lx05 = 0.015; 		// half Cell size in x (m) to set periodicity
double Ly05 = 0.015;            // half Cell size in y (m) to set periodicity

// double xC_fori[numero_fori] = {0.0}; // center of each aperture in x 
// double yC_fori[numero_fori] = {0.0}; // center of each aperture in y 



//mesh 
double mesh_x = 0.03;                                 //Mesh x_dimension [m]
double mesh_y = 0.03;                               //Mesh y_dimension [m]
double mesh_z = 0.554;                               //Mesh z_dimension [m]



// VOltages and Current
double V0=0.0;// grid 0 voltage in Volts (tipically PG)



const int nscan=1;// number of simulations 
double current_density_scan[nscan]={175}; // extracted currents, A/m2
double V1_scan[nscan]={5}; // grid 1 voltage (tipically EG)
double V2_scan[nscan]={115};// grid 2 voltage (tipically AG1/GG)
double V3_scan[nscan]={224};// grid 2 voltage (tipically AG1/GG)
double V4_scan[nscan]={334};// grid 2 voltage (tipically AG1/GG)
double V5_scan[nscan]={443};// grid 2 voltage (tipically AG1/GG)
double V6_scan[nscan]={552};// grid 2 voltage (tipically AG1/GG)


/*
const int nscan=5;// number of simulations 
double current_density_scan[nscan]={175, 175, 175, 175, 175}; // extracted currents, A/m2
double V1_scan[nscan]={5.2, 5.4, 5.6, 5.8, 6.0}; // grid 1 voltage (tipically EG)
double V2_scan[nscan]={115, 115, 115, 115, 115};// grid 2 voltage (tipically AG1/GG)
double V3_scan[nscan]={224, 224, 224, 224, 224};// grid 2 voltage (tipically AG1/GG)
double V4_scan[nscan]={334, 334, 334, 334, 334};// grid 2 voltage (tipically AG1/GG)
double V5_scan[nscan]={443, 443, 443, 443, 443};// grid 2 voltage (tipically AG1/GG)
double V6_scan[nscan]={552, 552, 552, 552, 552};// grid 2 voltage (tipically AG1/GG)
*/

/* grid voltages and current density used for the scans
const int nscan=24;// number of simulations 
double current_density_scan[nscan]={330,330,330,330,330,330,330,330,165,165,165,165,165,165,165,165,82,82,82,82,82,82,82,82}; // extracted currents, A/m2
double V1_scan[nscan]={8.8,9.0,9.2,9.4,9.6,9.8,10.0,10.2,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,3.5,3.55,3.6,3.65,3.7,3.75,3.8,3.85}; // grid 1 voltage (tipically EG)
double V2_scan[nscan]={209,209,209,209,209,209,209,209,127,127,127,127,127,127,127,127,78,78,78,78,78,78,78,78};// grid 2 voltage (tipically AG1/GG)
double V3_scan[nscan]={409,409,409,409,409,409,409,409,253,253,253,253,253,253,253,253,156,156,156,156,156,156,156,156};// grid 2 voltage (tipically AG1/GG)
double V4_scan[nscan]={609,609,609,609,609,609,609,609,379,379,379,379,379,379,379,379,234,234,234,234,234,234,234,234};// grid 2 voltage (tipically AG1/GG)
double V5_scan[nscan]={809,809,809,809,809,809,809,809,505,505,505,505,505,505,505,505,311,311,311,311,311,311,311,311};// grid 2 voltage (tipically AG1/GG)
double V6_scan[nscan]={1009,1009,1009,1009,1009,1009,1009,1009,631,631,631,631,631,631,631,631,390,390,390,390,390,390,390,390};// grid 2 voltage (tipically AG1/GG)
*/



/*const int nscan=8;// number of simulations 
double current_density_scan[nscan]={10,10,10,10,10,10,10,10}; // extracted currents, A/m2
double V1_scan[nscan]={700,750,800,850,900,1000,1050,1100}; // grid 1 voltage (tipically EG)

double V2_scan[nscan]={20000,20000,20000,20000,20000,20000,20000,20000,};// grid 2 voltage (tipically AG1/GG)
double V3_scan[nscan]={39000,39000,39000,39000,39000,39000,39000,39000};// grid 2 voltage (tipically AG1/GG)
double V4_scan[nscan]={58000,58000,58000,58000,58000,58000,58000,58000};// grid 2 voltage (tipically AG1/GG)
double V5_scan[nscan]={78000,78000,78000,78000,78000,78000,78000,78000};// grid 2 voltage (tipically AG1/GG)
double V6_scan[nscan]={97000,97000,97000,97000,97000,97000,97000,97000};// grid 2 voltage (tipically AG1/GG)

*/



/*
double V3_scan[nscan]={409000,409000,409000,409000,409000};// grid 4 voltage (tipically AG1/GG)
double V4_scan[nscan]={609000,609000,609000,609000,609000};// grid 5 voltage (tipically AG1/GG)
double V5_scan[nscan]={809000,809000,809000,809000,809000};// grid 6 voltage (tipically AG1/GG)
double V6_scan[nscan]={1009000,1009000,1009000,1009000,1009000};// grid 7 voltage (tipically AG1/GG)
*/



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
