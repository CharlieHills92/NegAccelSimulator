#ifndef _MITICA_GEOM2
#define _MITICA_GEOM2

//#include "../param/TEST_parameters.h"
#include "../globals.h"
#include "geom_function.h"
#include "func_solid.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

////////////////// Geometry of the electrodes  /////////////////////////

/// MITICA Accelerator, RFX 
// P. Veltri, IO,  01/2019.
// C. Poggi, IO, 02/2025

class MITICA_Accelerator {

    public:
    
	MITICA_Accelerator(void) {
        calc_grid_location();
    };
	MITICA_Accelerator(double ext_gap, double acc_gap) {
        set_extgridgap(ext_gap);
        set_accgridgap(acc_gap);
        calc_grid_location();
    };

    ~MITICA_Accelerator(){};

    void set_extgridgap(double ext_gap) {
        extgrid_gap = ext_gap;
    };
    void set_accgridgap(double acc_gap) {
        accgrid_gap = acc_gap;
    };

    void calc_grid_location(void) {

        double EGshift = PGz.back()+extgrid_gap;
        for(double& d : EGz) {d += EGshift;};

        double AG1shift = EGz.back()+accgrid_gap;
        for(double& d : AG1z) {d += AG1shift;};
        double AG2shift = AG1z.back()+accgrid_gap;
        for(double& d : AG2z) {d += AG2shift;};
        double AG3shift = AG2z.back()+accgrid_gap;
        for(double& d : AG3z) {d += AG3shift;};
        double AG4shift = AG3z.back()+accgrid_gap;
        for(double& d : AG4z) {d += AG4shift;};
        double AG5shift = AG4z.back()+accgrid_gap;
        for(double& d : AG5z) {d += AG5shift;};

    };

    // definition of PG grid
    bool PG_MITICA(double x, double y, double z) {
        double rc_start = 0.;double rc_end = 0.;
        return createGrid(x,y,z,PGz,PGr,rc_start,rc_end);
    };
    bool solid0_MITICA(double x, double y, double z) {return PG_MITICA(x,y,z);};

    // definition of EG grid
    bool EG_MITICA(double x, double y, double z) {
        double rc_start = 1.e-3;double rc_end = 1.7e-3;
        return createGrid(x,y,z,EGz,EGr,rc_start,rc_end);
    };
    bool solid1_MITICA(double x, double y, double z) {return EG_MITICA(x,y,z);};

    // definition of AG1 grid
    bool AG1_MITICA(double x, double y, double z) {
        double rc_start = 1.e-3;double rc_end = 2.e-3;
        return createGrid(x,y,z,AG1z,AG1r,rc_start,rc_end);
    };
    bool solid2_MITICA(double x, double y, double z) {return AG1_MITICA(x,y,z);};
    // definition of AG2 grid
    bool AG2_MITICA(double x, double y, double z) {
        double rc_start = 1.e-3;double rc_end = 2.e-3;
        return createGrid(x,y,z,AG2z,AG2r,rc_start,rc_end);
    };
    bool solid3_MITICA(double x, double y, double z) {return AG2_MITICA(x,y,z);};
    // definition of AG3 grid
    bool AG3_MITICA(double x, double y, double z) {
        double rc_start = 1.e-3;double rc_end = 2.e-3;
        return createGrid(x,y,z,AG3z,AG3r,rc_start,rc_end);
    };
    bool solid4_MITICA(double x, double y, double z) {return AG3_MITICA(x,y,z);};
    // definition of AG4 grid
    bool AG4_MITICA(double x, double y, double z) {
        double rc_start = 1.e-3;double rc_end = 2.e-3;
        return createGrid(x,y,z,AG4z,AG4r,rc_start,rc_end);
    };
    bool solid5_MITICA(double x, double y, double z) {return AG4_MITICA(x,y,z);};
    // definition of AG5 grid
    bool AG5_MITICA(double x, double y, double z){
        double rc_start = 1.e-3;double rc_end = 1.e-3;
        return createGrid(x,y,z,AG5z,AG5r,rc_start,rc_end);
    };
    bool solid6_MITICA(double x, double y, double z) {return AG5_MITICA(x,y,z);};

    // Create a FuncSolid object with a lambda capturing `this`
    FuncSolid create_PGsolid() {
        return FuncSolid([this](double x, double y, double z) {
            return this->PG_MITICA(x, y, z);
        });
    };
    FuncSolid create_EGsolid() {
        return FuncSolid([this](double x, double y, double z) {
            return this->EG_MITICA(x, y, z);
        });
    };

    FuncSolid create_AG1solid() {
        return FuncSolid([this](double x, double y, double z) {
            return this->AG1_MITICA(x, y, z);
        });
    };

    FuncSolid create_AG2solid() {
        return FuncSolid([this](double x, double y, double z) {
            return this->AG2_MITICA(x, y, z);
        });
    };


    FuncSolid create_AG3solid() {
        return FuncSolid([this](double x, double y, double z) {
            return this->AG3_MITICA(x, y, z);
        });
    };


    FuncSolid create_AG4solid() {
        return FuncSolid([this](double x, double y, double z) {
            return this->AG4_MITICA(x, y, z);
        });
    };


    FuncSolid create_AG5solid() {
        return FuncSolid([this](double x, double y, double z) {
            return this->AG5_MITICA(x, y, z);
        });
    };


    private:

    double accgrid_gap = 88e-3;
    double extgrid_gap = 6e-3;
    // PG coordinates
    std::vector<double> PGz_rel = {2.232986e-3,7e-3,7.4e-3,7.8e-3,9e-3};
    std::vector<double> PGr = {11e-3,7e-3,7e-3,7.5e-3,8e-3};
    // EG coordinates
    std::vector<double> EGz_rel = {0e-3,1e-3,6e-3,15.60411e-3,17e-3};
    std::vector<double> EGr = {7.5e-3,6.5e-3,6.5e-3,8.246201e-3,9.91878e-3};
    // AG1 coordinates
    std::vector<double> AG1z_rel = {0e-3,1e-3,15e-3,17e-3};
    std::vector<double> AG1r = {8e-3,7e-3,7e-3,9e-3};
    // AG2 coordinates
    std::vector<double> AG2z_rel = {0e-3,1e-3,15e-3,17e-3};
    std::vector<double> AG2r = {8e-3,7e-3,7e-3,9e-3};
    // AG3 coordinates
    std::vector<double> AG3z_rel = {0e-3,1e-3,15e-3,17e-3};
    std::vector<double> AG3r = {9e-3,8e-3,8e-3,10e-3};
    // AG4 coordinates
    std::vector<double> AG4z_rel = {0e-3,1e-3,15e-3,17e-3};
    std::vector<double> AG4r = {9e-3,8e-3,8e-3,10e-3};
    // AG5-GG coordinates
    std::vector<double> AG5z_rel = {0e-3,1e-3,16e-3,17e-3};
    std::vector<double> AG5r = {9e-3,8e-3,8e-3,9e-3};

    std::vector<double> PGz = PGz_rel;
    std::vector<double> EGz = EGz_rel;
    std::vector<double> AG1z = AG1z_rel;
    std::vector<double> AG2z = AG2z_rel;
    std::vector<double> AG3z = AG3z_rel;
    std::vector<double> AG4z = AG4z_rel;
    std::vector<double> AG5z = AG5z_rel;

};

////// End of Geometry 
////////////////////////////////////////////////////////////////////////////////

#endif


