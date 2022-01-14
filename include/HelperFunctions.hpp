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
 * File: HelperFunctions.hpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains all declarations for functions in the
 * HelperFunctions.cpp file. 
 *============================================================================*/

#ifndef HelperFunctions_hpp
#define HelperFunctions_hpp

#include <stdio.h>
#include "boost/multi_array.hpp"
using namespace std;
typedef boost::multi_array<double,3> multiarray;

//------------------------Struct Defs-------------------------------------------------------
struct sOptimalValues{
    double valueFn;
    double oc;
    double dNextOfficial;
    double vNextOfficial;
    double tNextOfficial;
    
    double valueFnParab;
    double ocParab;
};

struct sGSS{
    double valueFn;
    double oc;
    double dNext;
    double vNext;
    double tNext;
};

struct sOptimalTrajectory{
    vector<double> positionVals;
    vector<double> velocityVals;
    vector<double> timeVals;
    vector<double> optimalControlVals;
    vector<double> analyticalOCVals;
    
    //Vectors for Model 1 Tracing with two possible light lengths
    vector<double> s1Position;
    vector<double> s1Velocity;
    vector<double> s1Accel;
    vector<double> s1Time;
    double d0S1;
    double v0S1;
    
    vector<double> s2Position;
    vector<double> s2Velocity;
    vector<double> s2Accel;
    vector<double> s2Time;
    double d0S2;
    double v0S2;
    
    vector<double> s3Position;
    vector<double> s3Velocity;
    vector<double> s3Accel;
    vector<double> s3Time;
    double d0S3;
    double v0S3;
    
    vector<double> s4Position;
    vector<double> s4Velocity;
    vector<double> s4Accel;
    vector<double> s4Time;
    double d0S4;
    double v0S4;
    
    short int loopIdx;
};


/*==============================================================================
 * Loop over Controls
 *============================================================================*/
struct sOptimalValues controlLoop(multiarray *aU, multiarray *aOC, multiarray *aUYRG, multiarray *aOCYRG, multiarray *aUFirstPassYRG, multiarray *aOCFirstPassYRG, double aCurrentPos, double aCurrentVel, double aCurrentTime, short int aCurrentTimeIndex, double aTau, double aProbOfLightChange);

/*==============================================================================
* Golden Section Search (GSS)
*============================================================================*/
struct sGSS goldenSection (multiarray *aValueFunction, multiarray *aOC, multiarray *aUYRG, multiarray *aOCYRG, multiarray *aUYRGFirstPass, multiarray *aOCYRGFirstPass, short int aTimeIndex, double aLeftEndPt, double aRightEndPt, double aCurrentTime, double aDistCurrent, double aVelCurrent, double aTau);

/*==============================================================================
* Determination if value function is prescribed
*============================================================================*/
double prescribedValueFunction(multiarray *aU, multiarray *aOC, double aD, double aV, double aT, double aTau, double aControl);

/*==============================================================================
* Deterministic Running Cost Integral
*============================================================================*/
double runningCostIntegral(double aControl, double aCurrentTime, double aTau, double aGamma, double aC1, double aC2, double aC3);

/*==============================================================================
* Boundary cost along parabolic bd
*============================================================================*/
double parabBoundaryCost(multiarray *aU, multiarray *aOC, double aIntersectVel, double aIntersectT, double aTurningGreenT, double aGamma, double aC1, double aC2, double aC3, double aDeltaD, double aDeltaV, double aTargetPos);

/*==============================================================================
* Boundary cost along max accel line
*============================================================================*/
double maxAccelBoundaryCost(multiarray *aU, multiarray *aOC, double aIntersectVel, double aIntersectT, double aTurningRedT, double aGamma, double aC1, double aC2, double aC3, double aDeltaD, double aDeltaV, double aTargetPos);

/*==============================================================================
* Boundary cost along target at d = df
*============================================================================*/
double targetBoundaryCost(double aV, double aT);

/*==============================================================================
* Terminal cost
*============================================================================*/
double terminalCost (multiarray *aU, multiarray *aOC, multiarray *aUParab, multiarray *aOCParab, double aD, double aV, double aT, double aCurrentPos, double aCurrentVel, double aCurrentTime, double aControl, double aMu, double aTargetPos);

/*==============================================================================
 * Trajectory tracer
 *============================================================================*/
struct sOptimalTrajectory trajectoryTracer(multiarray *aValueFunction, multiarray *aOptimalControl, multiarray *aValueFunctionYRG, multiarray *aOptimalControlsYRG, multiarray *aValueFunctionYRGFirstPass, multiarray *aOptimalControlsYRGFirstPass, double aPosStart, double aVelStart, double aTimeStart, double aTau);

/*==============================================================================
 * All trajectories for Model 1, Two Lengths
 *============================================================================*/
void model1TwoLengthsTrajectories(multiarray *aValueFunctionUncertainGreen, multiarray *aOCUncertainGreen, multiarray *aYRGValueFunction, multiarray *aYRGOC, multiarray *aYRGFirstPassValueFunction, multiarray *aYRGFirstPassOC, multiarray *aStationaryValueFunction, multiarray *aStationaryOC, double aPosStart, double aVelStart, double aTimeStart, double aTau);

/*==============================================================================
 * All trajectories for Model 1, Three Lengths
 *============================================================================*/
void model1ThreeLengthsTrajectories(multiarray *aValueFunctionUncertainGreen, multiarray *aOCUncertainGreen, multiarray *aYRGValueFunction, multiarray *aYRGOC, multiarray *aYRGFirstPassValueFunction, multiarray *aYRGFirstPassOC, multiarray *aStationaryValueFunction, multiarray *aStationaryOC, double aPosStart, double aVelStart, double aTimeStart, double aTau);

/*==============================================================================
 * All trajectories for Model 1, Four Lengths
 *============================================================================*/
void model1FourLengthsTrajectories(multiarray *aValueFunctionUncertainGreen, multiarray *aOCUncertainGreen, multiarray *aYRGValueFunction, multiarray *aYRGOC, multiarray *aYRGFirstPassValueFunction, multiarray *aYRGFirstPassOC, multiarray *aStationaryValueFunction, multiarray *aStationaryOC, double aPosStart, double aVelStart, double aTimeStart, double aTau);

/*==============================================================================
 * Pareto Grid
 *============================================================================*/
void paretoGrid(short int aNumPointsADim, short int aNumPointsBDim, double aC1, double aC2, double aC3, double aTrajStartPos, double aTrajStartVel, double aTrajStartTime);

/*==============================================================================
 * J3 Constrained Pareto Front Generation
 *============================================================================*/
void j3ConstrainedParetoFront(short int aNumPointsADim, short int aNumPointsBDim, double aJ3ConstantValue, double aTrajStartPos, double aTrajStartVel, double aTrajStartTime);

/*==============================================================================
 * Mollifying Terminal Condition
 *============================================================================*/
void mollifyTerminalCondition (multiarray *aInitialValues, multiarray *aMollifiedSoln, short int aNumberOfTimesteps, double aDiffusionCoef);

/*==============================================================================
 * Infinity Norm
 *============================================================================*/
double infinityNorm(multiarray *aU, short int aPosIdx);

/*==============================================================================
 * Write to file
 *============================================================================*/
void writeToFile (multiarray *aU, multiarray *aControls);

/*=======================================================================================================================
 * Write to file - row major order
 *======================================================================================================================*/
void rowMajorWriteToFile(multiarray *aValueFunction, multiarray *aOptimalControls);

/*==============================================================================
 * Write to file
 *============================================================================*/
void readFromFile(multiarray *aValueFunctionGRG, multiarray *aOptimalControlsGRG);

/*==============================================================================
 * Write to file - trajectory
 *============================================================================*/
void trajWriteToFile (vector<double> aPosVec, vector<double> aVelVec, vector<double> aTimeVec, vector<double> aOCVals, short int aLoopEndIdx);

/*==============================================================================
 * Write to file - pareto
 *============================================================================*/
void paretoWriteToFile (multiarray *aJ1Array, multiarray *aJ2Array, multiarray *aJ3Array, vector<double> aDiscomfortCosts, vector<double> aFuelCosts, vector<double> aTimeCosts, vector<double> aC1Vals, vector<double> aC2Vals, vector<double> aC3Vals, short int aNumAPoints, short int aNumBPoints);

/*==============================================================================
 * read from file - pareto
 *============================================================================*/
void paretoReadFromFile(multiarray *aJ1Vals, multiarray *aJ2Vals, multiarray *aJ3Vals, short int aNumPointsADim, short int aNumPointsBDim);

/*==============================================================================
 * Write to file - debug arrays
 *============================================================================*/
void debugWriteToFile (vector<double> aControlVec, vector<double> aVFVec);

#endif /* HelperFunctions_hpp */
