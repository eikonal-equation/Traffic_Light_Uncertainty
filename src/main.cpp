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

//=======================================Libraries=======================================
#include <iostream>
#include <fstream>
#include <cmath>

//==================================Project-specific Header files=========================
#include "boost/multi_array.hpp"
#include "variables.h"
#include "InitFunctions.hpp"
#include "ProblemFunctions.hpp"
#include "HelperFunctions.hpp"

/*==============================================================================
 * File: main.cpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file initializes global variables and
 * executes runOptimalDrivingExamples() from the command line.
 *============================================================================*/

/*==============================================================================
 * CONSTANTS AND PARAMETER DEFINITIONS -
 *
 * Initialize the constants and parameters used
 * throughout the simulations.
 *
 * Note: Many of these values are set in initGrid().
 *============================================================================*/
bool gArtificial = false; //Toggle between artifical and realistic parameters

const long gInfty = 10000000; // Approximation of infinity
double gAccelInfty; //Acceleration value when u=infinity
double gVMax; // max velocity
double gDMax; // max distance

double gDTarget; //location of target
double gTerminalT; //just a placeholder for now

//==========================Controls and Maximum Red Light Duration===========================================
double gAlpha; //artificial magnitude of the maximum deceleration
double gBeta; //artificial magnitude of max accel
double gRedLightLength = 1; //Artificial time length of red light. Not currently in use

//Setting considered controls based on GSS or grid
bool gGoldenSection = true;
int gControlGridDim = 21; //number of possible control values to consider between -alpha, 0
int gNumControls = gControlGridDim; //set number of controls based on whether GSS used or not
vector<double> gControlValues; //vector of control values
double gControlSpacing; //spacing between control values on the grid

//==========================Parameters for grid creation=========================================================
double dtDenom; // max velocity and max accel magnitude for restrictions
int gN; // starting number of distance and velocity intervals
double gH; //Position and velocity spacing assuming that gDMax = gVMax //debugging experiment

double gDT; //set DT based on timestep restriction
double gDPos; //init delta D
double gDVel;// init delta V
double gTau; //default starting tau
int gNt; //number of time intervals
int gDNum; //number of distance intervals
int gVNum; //velocity intervals

//Counters to keep track of the i,j,k indices in time dependent solver loop. Must be init here to 0
int gTimeIndex = 0; //k
int gPosIndex = 0; //i
int gVelIndex = 0; //j

//Grid refinement
bool gGridRef; //grid ref switch
int gGridRefFactor; //grid ref factor

//=================================================Parameters for light length information======================================

double gTG; //red phase end time
double gTR; //red phase start time
short int gTGIdx = short(int(round(gTG / gDT))); //turning green index plceholder
short int gTimeEndIdx = 0; //Minimum index for the loop over time. Default set to 0.

short int gNumStops = 2; //number of possible light lengths - officialy set in initGrid()
vector<double> gLightTimes; //vector of turning yellow times
vector<double> gLightProbs; //vector of probabilities associated with each turning yellow time

double gTY1; //first turning yellow time
double gTY2; //second turning yellow time
double gTY3; //third turning yellow time
double gTY4; //fourth turning yellow time

double gTYMax; //maximum turning yellow time
double gTYP1; //probability of TY1
double gTYP2; //probability of TY2
double gTYP3; //probability of TY3
double gTYP4; //probability of TY4

double gYellowDuration; //length of yellow light (seconds)
double gRedDuration; //length of red light (seconds)

//==============================Objective paramaters======================================================
double gGamma = 1;
double gMu = 1;

double gC1 = 1.0/3.0; //fuel cost objective weight
double gC2 = 1.0/3.0; //discomfort cost objective weight
double gC3 = 1.0/3.0; //time to target objective weight

//------------------------------------------Debugging-----------------------------------
bool gDebugControlChatter = false;

vector<double> gValueFunctionValsAtTestPoint;
vector<double> gTerminalCondValsAtTestPoint;
vector<double> gUNextValsAtTestPoint;
vector<double> gTimeVals;
vector<double> gProbVals;
vector<double> gTestParabXPts;
vector<double> gTestParabYPts;
vector<double> gParabTestXNodes; //holds dLeft, dPoint, dMiddle, and dRight
vector<double> gParabTestYNodes; //holds Q1, Q2, Q0, and vF at P
int gNumParabPts = 100;

bool gMollifyTerm = false;
bool gHorizDiffusionProbOnly = true;

bool gTrajTracerDebug = false;
double gVelocityTimeDebugCounter = 0;

//======================================Sampling for Larger Grid Sizes - NOT CURRENTLY USED==================================
short int gSamplingTNumThresh = 10000; //minimum number of time intervals needed to start sampling procedure instead of writing all to file
short int gSliceSampleSpacing = 10; //sample every 10 intervals

//========================================PROBLEM SELECTION SWITCHES==========================================
bool gInStationarySolve; //indicates if solving stationary problem
bool gSolvingGreenOnly; //indicates if in Example 1
bool gInGS; //indicates if in Gauss-Seidel routine
bool gModifiedGS; //indicates if in Modifiec Gauss-Seidel routine
bool gInGJ; //indicates if in Gauss-Jordan routine

bool gGreenRedGreenProb; //indicates if in Example 2 or 3
bool gInTSolve; //indicates if solving for time to target PDE (Not currently used)
bool gInRGSolve; //Indicates if solving RG ONLY problem

bool gUncertainGreenProb; //indicates if in Example 4 or 5
bool gInUncertainGreenPhase; //indicates if we are solving the PDE in the uncertain initial green phase
bool gInTerminalPenalty; //indicates if we are computing the terminal cost in uncertain initial green phase

//=========================================UNCERTAINTY SWITCHES / INIT========================================
bool gModel1; //indicates if in model 1 uncertainty
bool gModel2; //indicates if in model 2 uncertainty (not currently used)
bool gPedestrianModel; //indicates if in pedestrian actuated uncertainty (not currently used)

double gArrivalRate = (5.0/60.0); //pedestrian arrival rate for pedestrian actuated uncertainty model (not currently used)

bool gInEqualCEqualP2lExample; //indicates if in example where c's equal and probs equal for Example 4
bool gInC2Priority2lExample; //indicates if in example where c2 has highest priority for Example 4
//========================================INTERPOLATION SWITCHES==============================================
bool gHasBoundaryOn; //indicates if piecewise boundary in effect
bool gFirstPass; //indicates if in first sweep of YRG problem
bool gInMaxAccelInterp; //indicates if in interpolation function along max accel line
bool gInParabInterp;//indicates if in parabolic interpolation function

//=========================================TRAJECTORY TRACER SWITCHES=========================================
bool gInTracer; //indicates if in trajectory tracer
bool gControlInterp; //indicates if using control interpolation in tracer
bool gStationaryTraj; //indicates if tracing trajectory for the stationary problem

bool gModel1Traj; //indicates if model 1 uncertainty used and model 1 trajectory is being traced
bool gModel1TrajDetPortion; //indicates if tracing deterministic portion of a model 1 trajectory
bool gFirstTraj; //indicates if tracing trajectory associated with turning yellow at T_1 (Examples 4 and 5)
bool gSecondTraj; //indicates if tracing trajectory associated with turning yellow at T_2 (Examples 4 and 5)
bool gThirdTraj; //indicates if tracing trajectory associated with turning yellow at T_3 (Example 5)
bool gFourthTraj; //indicates if tracing trajectory associated with turning yellow at T_4 (Not used in any examples)

//======================================PARETO SWITCHES=======================================================
bool gPareto; //indicates if we are finding Pareto front
bool gSelectedTraj; //indicates trajectory tracing in Pareto framework
bool gParetoCostComputation; //indicates if we are computing cost along trajectory in Pareto
bool gStationaryPareto; //indicates if finding Pareto front during indefinite green (Example 1)
bool gYRGPareto; //indicates if finding pareto front during YRG (not currently used)
bool gUncertainPareto; //indicates if finding Pareto front during uncertain green (not currently used)
bool gInConstrainedPareto; //indicates if finding constrained Pareto front (Example 1)
bool gJ3ConstrainedCurve1; //indicates if finding j3-constrained Pareto for first time to target tolerance (Example 1)
bool gJ3ConstrainedCurve2; //indicates if finding j3-constrained Pareto for second time to target tolerance (Example 1)
bool gJ3ConstrainedCurve3; //indicates if finding j3-constrained Pareto for third time to target tolerance (Example 1)

//Costs along trajectory in Pareto computation. Must be init to 0 and reset after each trajectory trace in Pareto
double gDiscomfortCost = 0; //discomfort cost along trajectory
double gFuelCost = 0; //fuel cost along trajectory
double gTimeCost = 0; //time to target
double gTotalCost = 0; //total cost
//===================================EXAMPLE SWITCHES=========================================================
bool gInExample1TrajectoryComputation = false; //indicates if we are computing the sample trajectory in Example 1

//==============================MAIN IMPLEMENTATION============================================================
int main(int argc, const char * argv[]) {
    /*
     Command line argument key:

     argv[1] = example number.

     */

    setFrameworkSettings();

    int problemSelector = stoi(argv[1]);

    std::cout << "Example selected: " << argv[1] << "\n";

    runOptimalDrivingExamples(problemSelector);

    std::cout << "Computation complete. Results saved to file(s). \n";

    return 0;
}
