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
 * File: InitFunctions.hpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains the function declarations of functions utilized
 * in initializing the grids, parameters, and variables used throughout Optimal Driving.
 * Associated functions are defined and described in InitFunctions.cpp
 *============================================================================*/
#ifndef InitFunctions_hpp
#define InitFunctions_hpp

#include <stdio.h>

//using namespace std;
typedef boost::multi_array<double,3> multiarray;

/*==============================================================================
 * Set booleans
 *============================================================================*/
void setFrameworkSettings();

/*==============================================================================
 * Init grid ob extension
 *============================================================================*/
void initGrid();

/*==============================================================================
 * Control Array Init
 *============================================================================*/
void initControls(int aControlArrayDim);

/*==============================================================================
 * Array initialization
 *============================================================================*/
//void initArray(multiarray *aValueFnArray, multiarray *aOptimalControlArray, const double aOptControl);

/*==============================================================================
 * Grid refinement
 *============================================================================*/
void gridRef (double aDeltaD, double aDeltaV, double aDeltaT, double aRefinementFactor);

#endif /* InitFunctions_hpp */

