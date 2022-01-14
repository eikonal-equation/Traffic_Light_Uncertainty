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
 * File: StationarySolver.hpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains the declarations for the functions which solve
 * the stationary HJB equation using Gauss-Jacobi and Gauss-Seidel iterations.
 *
 * The corresponding functions are described in StationarySolver.cpp
 *============================================================================*/
#ifndef StationarySolver_hpp
#define StationarySolver_hpp

#include <stdio.h>
#include "boost/multi_array.hpp"

using namespace std;
typedef boost::multi_array<double,3> multiarray;

/*==============================================================================
 * Stationary Green Light Solver - Gauss Jacobi
 *============================================================================*/
void stationaryGJSolver(multiarray *aU, multiarray *aOC);

/*==============================================================================
 * Stationary Green Light Solver - Gauss Seidel
 *============================================================================*/
void stationaryGSSolver(multiarray *aU, multiarray *aOC);

#endif /* StationarySolver_hpp */

