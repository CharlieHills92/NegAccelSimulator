#include "MITICA_geom2.h"

MITICA_Accelerator::~MITICA_Accelerator()
{
}

// PG grid definition
bool MITICA_Accelerator::PG_MITICA(double x, double y, double z) {
    double rc_start = 0.;double rc_end = 0.;
    return createGrid(x,y,z,PGz,PGr,rc_start,rc_end);
}

// EG grid definition
bool MITICA_Accelerator::EG_MITICA(double x, double y, double z) {
    double rc_start = 1.e-3;double rc_end = 1.7e-3;
    return createGrid(x,y,z,EGz,EGr,rc_start,rc_end);
}

// AG1 grid definition
bool MITICA_Accelerator::AG1_MITICA(double x, double y, double z) {
    double rc_start = 1.e-3;double rc_end = 2.e-3;
    return createGrid(x,y,z,AG1z,AG1r,rc_start,rc_end);
}

// AG1 grid definition
bool MITICA_Accelerator::AG2_MITICA(double x, double y, double z) {
    double rc_start = 1.e-3;double rc_end = 2.e-3;
    return createGrid(x,y,z,AG2z,AG2r,rc_start,rc_end);
}

// AG1 grid definition
bool MITICA_Accelerator::AG3_MITICA(double x, double y, double z) {
    double rc_start = 1.e-3;double rc_end = 2.e-3;
    return createGrid(x,y,z,AG3z,AG3r,rc_start,rc_end);
}

// AG1 grid definition
bool MITICA_Accelerator::AG4_MITICA(double x, double y, double z) {
    double rc_start = 1.e-3;double rc_end = 2.e-3;
    return createGrid(x,y,z,AG4z,AG4r,rc_start,rc_end);
}

// AG1 grid definition
bool MITICA_Accelerator::AG5_MITICA(double x, double y, double z) {
    double rc_start = 1.e-3;double rc_end = 1.e-3;
    return createGrid(x,y,z,AG5z,AG5r,rc_start,rc_end);
}

