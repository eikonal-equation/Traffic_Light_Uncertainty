/*=====================================================================================
 * Copyright (C) 2021 Mallory E. Gaspard
 *
 * This program is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see http://www.gnu.org/licenses/.
 *===================================================================================*/

/*==============================================================================
 * File: ProbabilityFunctions.hpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains function declarations for the functions
 * described in ProbabilityFunctions.cpp.
 *============================================================================*/
#ifndef ProbabilityFunctions_hpp
#define ProbabilityFunctions_hpp

#include <stdio.h>

/*==============================================================================
 * Conditional probability determination
 *============================================================================*/
double conditionalProb(double aT);

/*==============================================================================
 * Prob of light change for discrete turning yellow times
 *============================================================================*/
double discreteProb(double aT);

/*==============================================================================
 * Prob of light change for uniform light length problem
 *============================================================================*/
double uniformProb(double aT);

#endif /* ProbabilityFunctions_hpp */
