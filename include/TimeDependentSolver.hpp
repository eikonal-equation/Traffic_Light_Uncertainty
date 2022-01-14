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
 * File: TimeDependentSolver.hpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains the declarations of the function responsible for executing
 * the time-dependent HJB solve during the red phase, yellow phase, and the
 * initial uncertain green phase and the function responsible for executing the second-pass
 * sweep during the yellow phase in the YRG problem to account for any discontinuity
 * present. Functions themselves are in the associated cpp file.
 *============================================================================*/

#ifndef TimeDependentSolver_hpp
#define TimeDependentSolver_hpp

#include <stdio.h>
#include "boost/multi_array.hpp"
using namespace std;
typedef boost::multi_array<double,3> multiarray;

/*==============================================================================
* Yellow-Red Solver
*============================================================================*/
void timeDepSolver(multiarray *aU, multiarray *aOC, multiarray *aUYRG, multiarray *aOCYRG, multiarray *aUYRGFirstPass, multiarray *aOCYRGFirstPass, multiarray *aStationarySoln, multiarray *aStationaryOC);
/*==============================================================================
 * Value function corrector for [TY, TR]
 *============================================================================*/
void valueFunctionYellowPhaseResolve(multiarray *aValueFunction, multiarray *aOptimalControls, multiarray *aValueFunctionFirstPass, multiarray *aOptimalControlsFirstPass, multiarray *aStationarySoln, multiarray *aStationaryOC);

#endif /* TimeDependentSolver_hpp */

