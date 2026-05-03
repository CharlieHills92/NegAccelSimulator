/*
 * StrippingUtils.cpp
 *
 *  Created on: Aug 07, 2025
 *      Author: GitHub Copilot
 */

#include "StrippingUtils.h"

std::string getDensityProfileFilename(uint accelerator_idx) {
    switch (accelerator_idx) {
        case 1: // SPIDER - use MITICA as fallback
            return "densprofiles/MITICA_dens.txt";
        case 2: // MITICA
            return "densprofiles/MITICA_dens.txt";
        case 3: // MTF
            return "densprofiles/MTF_dens.txt";
        default:
            return "densprofiles/MITICA_dens.txt"; // Default fallback
    }
}
