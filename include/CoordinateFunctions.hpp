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
 * File: CoordinateFunctions.hpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains the declarations of the tau refinement function,
 * as well as the declarations of the functions responsible for determining whether
 * a (d,v) point belongs to the allowable region at time t, where a (d,v) point
 * lies in relation to the pw boundary at time t, the location of the maximum
 * acceleration line d_{\beta} at time t, and the threshold velocity at time t.
 *
 * Details of all of these functions are found in CoordinateFunctions.cpp.
 *
 *============================================================================*/
#ifndef CoordinateFunctions_hpp
#define CoordinateFunctions_hpp

#include <stdio.h>

struct sNextCoordinates{
    double nextPos;
    double nextVel;
    double nextTime;
    double tauValue;
};

/*==============================================================================
* Determine dNext, vNext, tNext and tau refined
*============================================================================*/
struct sNextCoordinates tauRefinement(double aD, double aV, double aT, double aControl, double aTau);

/*==============================================================================
 * Determining if (d,v) point in allowable region at time t
 *============================================================================*/
bool isAllowed(double aV, double aD, double aT);

/*==============================================================================
 * Determining if d,v in allowable region during the yellow and red phases
 *============================================================================*/
bool isAllowedRedPhase(double aV, double aD, double aT);

/*=====================================================================================
 * Determining if starting d,v is in the legal region above max accel line
 *=======================================================================================*/
bool isInAccelRegion(double aV, double aD, double aT);

/*=======================================================================================
 * Location of (d,v) point in relation to PARABOLIC boundary d_{\alpha}
 *=====================================================================================*/
bool isInParabola(double aV, double aD, double aT);

/*==============================================================================
 * Is on considered portion of the parabolic bd d_{\alpha}
 *============================================================================*/
bool isOnParabPiece (double aD, double aV, double aT);

/*==============================================================================
 * Is on considered portion of the max accel bd d_{\beta}
 *============================================================================*/
bool isOnMaxAccelBD (double aD, double aV, double aT);

/*========================================================================================
 * Determining if starting d,v is subject to the max accel line or the extension
 *========================================================================================*/
bool isSubjectToMaxAccelExtension(double aV, double aD, double aT);

/*==============================================================================
 * Determine max accel line position for a given v
 *============================================================================*/
double maxAccelLine (double aV, double aT);

/*==================================================================================
 * Determining velocity at intersection between max accel line and parabolic bd
 *==================================================================================*/
double threshVelocityMaxAccel(double aD, double aV, double aT);

#endif /* CoordinateFunctions_hpp */

