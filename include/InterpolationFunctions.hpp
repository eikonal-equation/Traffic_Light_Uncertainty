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
 * File: InterpolationFunctions.hpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains all of the function declarations for the functions
 * defined in InterpolationFunctions.cpp
 *============================================================================*/
#ifndef InterpolationFunctions_hpp
#define InterpolationFunctions_hpp

#include <stdio.h>
#include "boost/multi_array.hpp"
using namespace std;
typedef boost::multi_array<double,3> multiarray;

struct sOptimalValuesControlInterp{
    double valueFn;
    double oc;
    double dNextOfficial;
    double vNextOfficial;
    double tNextOfficial;
};

/*==============================================================================
 * Bilinear Interpolation scheme
 *============================================================================*/
double bilinearInterp(multiarray *aU, multiarray * aOC, short int aTimeIndex, double aPos, double aVel, double aTime, double aControl);

/*======================================================================================
 * Cut cell bilinear Interpolation scheme for problems with max-accel bd
 *=====================================================================================*/
double maxAccelCutCellBilinearInterp (multiarray *aU, multiarray * aOC, short int aTimeIndex, double aNextPos, double aNextVel, double aNextTime, double aCurrentPos, double aCurrentVel, double aCurrentTime, double aControl);

/*==============================================================================
 * Cut cell bilinear Interpolation scheme for parabolic portion
 *============================================================================*/
double parabCutCellBilinearInterp (multiarray *aU, multiarray *aOC, short int aTimeIndex, double aNextPos, double aNextVel, double aNextTime, double aCurrentPos, double aCurrentVel, double aCurrentTime, double aControl);

/*============================================================================================================================
* Discontinuity Line Interp
*==========================================================================================================================*/
double discontinuityLineInterp(multiarray *aU, multiarray *aOC, short int aTimeIndex, double aNextPos, double aNextVel, double aNextTime, double aCurrentPos, double aCurrentVel, double aCurrentTime, double aControl);

/*==============================================================================
 * Entire interpolation routine
 *============================================================================*/
double generalInterp(multiarray *aU, multiarray * aOC, multiarray *aUParab, multiarray *aOCParab, short int aTimeIndex, double aNextPos, double aNextVel, double aNextTime, double aCurrentPos, double aCurrentVel, double aCurrentTime, double aControl);

/*==============================================================================
* Interpolation in d,v,t
*============================================================================*/
double dvtInterpolation(multiarray *aU, multiarray *aOC, multiarray *aUParab, multiarray *aOCParab, double aV, double aD, double aDistCurrent, double aVelCurrent, double aCurrentTime, double aTimeUpdated, double aControlValue);

/*==============================================================================
 * Interpolating between control values
 *============================================================================*/
struct sOptimalValuesControlInterp controlInterp(multiarray *aU, multiarray *aOC, multiarray *aUFirstPass, multiarray *aOCFirstPass, double aCurrentPos, double aCurrentVel, double aCurrentTime, short int aCurrentTimeIndex, double aTau, double aProbOfLightChange);

#endif /* InterpolationFunctions_hpp */

