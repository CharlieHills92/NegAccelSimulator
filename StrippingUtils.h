/*
 * StrippingUtils.h
 *
 *  Created on: Aug 07, 2025
 *      Author: GitHub Copilot
 */

#ifndef STRIPPINGUTILS_H_
#define STRIPPINGUTILS_H_

#include <string>

/**
 * @brief Get the appropriate density profile filename based on accelerator type
 * @param accelerator_idx Accelerator index (1=SPIDER, 2=MITICA, 3=MTF)
 * @return Path to the density profile file
 */
std::string getDensityProfileFilename(uint accelerator_idx);

#endif /* STRIPPINGUTILS_H_ */
