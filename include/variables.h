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


#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>

//Header Files
#include "boost/multi_array.hpp"

/*==============================================================================
 * File: variables.h
 *
 * Author: Mallory Gaspard
 *
 * Description: This file holds the definition of global variables used
 * throughout the Optimal Driving framework. See main.cpp for a description
 * of each of their functions.
 *============================================================================*/

#ifndef variables_h
#define variables_h

using namespace std;
typedef boost::multi_array<double,3> multiarray;

/*-----------------------------------------------------------------------------
 * Global Parameters:
 *----------------------------------------------------------------------------*/

/*-------------------------------------------Physical Params---------------------------------------------*/
extern const long gInfty; // Approximation of infinity
extern double gAccelInfty; //Acceleration value when u=infinity

extern double gVMax; // max velocity
extern double gDMax; // max distance
extern double gAlpha; //magnitude of the maximum decceleration
extern double gBeta; //max accel
extern double gRedLightLength; //Time length of red light

extern double gC1;
extern double gC2;
extern double gC3;
extern double gGamma; //fuel consumption parameter
extern double gMu; //terminal cost param in OB Extension problem
extern double gDTarget; //target distance
extern double gTerminalT; //terminal time
extern double gDPos; //pos spacing
extern double gDVel; //vel spacing

/*--------------------------------------------------Red Phase Parameters--------------------------------------------*/
extern double gTG; //ending time of red phase
extern double gTR; //starting time of red phase
extern short int gTGIdx;
extern short int gTimeEndIdx;
extern bool gHasBoundaryOn;
extern bool gGreenRedGreenProb;
extern bool gUncertainGreenProb;
extern bool gGridRef;
extern bool gQuadInterp;
/* -------------------------------------------Uncertain Green Phase--------------------------------------------*/
extern double gTY1;
extern double gTY2;
extern double gTY3;
extern double gTY4;
extern double gTYMax;
extern double gTYP1;
extern double gTYP2;
extern double gTYP3;
extern double gTYP4;
extern double gYellowDuration;
extern double gRedDuration;
extern bool gInUncertainGreenPhase;
extern bool gSolvingGreenOnly;
/*-----------------------------------------------------Controls---------------------------------------------------------------------------------------*/
extern int gNumControls; //number of control values considered
extern vector<double> gControlValues;
extern bool gGoldenSectionOn;
extern int gControlGridDim;
extern bool gParabBDIsExtended;
/*-----------------------------------------------Parameters for grid creation--------------------------------------------------------*/
extern bool gArtificial; //toggle artificial and realistic params
extern double gDT; //time-step size
extern double gH; // grid spacing
extern int gN; //number of spatial slices
extern int gNt; //number of timeslices
extern int gDNum; //number of distance slices
extern int gVNum;
extern double gControlSpacing; //spacing between gridded control values
extern double gTau; //SL step
extern int gTimeIndex;
extern int gPosIndex;
extern int gVelIndex;
extern int gGridRefFactor;
extern short int gIterations; //number of times to refine grid in convergence loop
extern bool gInStationarySolve;
extern bool gInGS;
extern bool gInGJ;
extern bool gModifiedGS;
extern bool gInTSolve;
extern bool gInRGSolve;

extern bool gFirstPass; //for the [0, TR] correction
extern bool gDiscontLineInterpUsed;
extern bool gInMaxAccelInterp;
extern bool gInParabInterp;
extern bool gUpBiased;
extern bool gDownBiased;
extern bool gRightBiased;
extern bool gLeftBiased;
extern bool gInTracer; //switch to indicate we're in the trajectory tracer
extern bool gControlInterp;
extern bool gGoldenSection; //turn on / off gss for controls
extern short int gNumStops; //number of possible light lengths (discrete)
extern vector<double> gLightProbs; //array to hold light length probabilities
extern vector<double> gLightTimes; //array to hold possible light length times

extern bool gPareto;
extern bool gStationaryPareto;
extern bool gYRGPareto;
extern bool gUncertainPareto;
extern bool gInConstrainedPareto;
extern bool gJ3ConstrainedCurve1;
extern bool gJ3ConstrainedCurve2;
extern bool gJ3ConstrainedCurve3;

extern bool gSelectedTraj;
extern bool gParetoCostComputation;
extern double gDiscomfortCost;
extern double gFuelCost;
extern double gTimeCost;
extern double gTotalCost; 
/*--------------------------------------------Probability / Model Info --------------------------------------------------------------------*/
extern bool gModel0;
extern bool gModel1;
extern bool gModel2;

extern bool gPedestrianModel;
extern double gArrivalRate;

extern bool gInEqualCEqualP2lExample;
extern bool gInC2Priority2lExample;
extern bool gInTerminalPenalty;
extern bool gModel1Traj;
extern bool gFirstTraj;
extern bool gSecondTraj;
extern bool gThirdTraj;
extern bool gFourthTraj;
extern bool gModel1TrajDetPortion;
extern bool gStationaryTraj;

extern bool gMollifyTerm;
extern bool gHorizDiffusionProbOnly;
/*-----------------------------------------------Example switches---------------------------------------------------------------------------*/
extern bool gInExample1TrajectoryComputation;
/*--------------------------------------------------------------Sampling for Larger Grids---------------------------------------------------*/
extern short int gSamplingTNumThresh; //minimum number of timeslices needed to start sampling procedure instead of writing all to file
extern short int gSliceSampleSpacing; //sample every so many slices


extern bool gTrajTracerDebug;
extern bool gDebugCutCells;
extern double gVelocityTimeDebugCounter;
extern bool gDebugControlChatter;
extern double gNegativeQuadCounter;
extern int gNegativeControlCounterGreenOnly;
extern bool gBilinearInterpUsed;
extern double gOriginalInterpVal;
extern int gNumBLUsedAndControlNegative;
extern bool gBilinearInterpUsedForOptimal;
extern bool gBilinearInterpUsedAtLeastOnce;
extern int gBLInterpUsedAtLeastOnceCounter;
extern bool gInCutCellDet;

extern vector<double> gValueFunctionValsAtTestPoint;
extern vector<double> gTerminalCondValsAtTestPoint;
extern vector<double> gUNextValsAtTestPoint;
extern vector<double> gTimeVals;
extern vector<double> gProbVals;
extern vector<double> gTestParabXPts;
extern vector<double> gTestParabYPts;
extern vector<double> gParabTestXNodes; //holds dLeft, dPoint, dMiddle, and dRight
extern vector<double> gParabTestYNodes; //holds Q1, Q2, Q0, and vF at P
extern int gNumParabPts;
extern bool gSwitchToBLInterp;
extern bool gInBilinearInterpRerun;
extern int gFeasibilityCondFailed;

#endif /* variables_h */
