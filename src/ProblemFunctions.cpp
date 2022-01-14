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
 * File: ProblemFunctions.cpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains all of the functions responsible for
 * executing each type of problem described in the Optimal Driving manuscript
 * as well as the four examples presented in Section 5 of the manuscript.
 *============================================================================*/

#include "ProblemFunctions.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include "variables.h"
#include "InitFunctions.hpp"
#include "StationarySolver.hpp"
#include "boost/multi_array.hpp"
#include "HelperFunctions.hpp"
#include "TimeDependentSolver.hpp"

using namespace std;

typedef boost::multi_array<double, 3> multiarray;

/*==============================================================================
 * Green Light Only Problem
 *============================================================================*/
void greenLightOnlyProblem(double aTrajStartPos, double aTrajStartVel, double aTrajStartTime){
    /*
     Function: greenLightOnlyProblem

     Arguments: starting (d,v,t) for optimal trajectory tracing.

     Purpose: Solves the stationary green light only problem using Modified Gauss Seidel.
     Traces optimal trajectory to target from supplied starting point.
     */

    gSolvingGreenOnly = true;
    gInStationarySolve = true;
    gGreenRedGreenProb = false;

    initGrid();

    short int tPosMaxIdx = (gDNum / 2) + 1;

    //Step 1: Set up d x v matrices to hold solution
    multiarray valueFunctionStationary(boost::extents[gDNum+1][gVNum+1][2]); //array for values of value function at each point in domain
    multiarray optimalControlsStationary(boost::extents[gDNum+1][gVNum+1][2]); //array for acceleration at each point in domain

    //Init array
    for(short int i = 0; i < gDNum + 1; i++) {     // distance loop
        for(short int j = 0; j < gVNum + 1; j++) { // velocity loop
            for(short int k = 0; k < 2; k++) {   // time loop

                double d = i * gDPos + gDTarget; //shifted for target at -1

                double distFromTarget = d - gDTarget; //target somewhere to left of d = 0
                double termCost = gMu * pow(distFromTarget, 2); //terminal cost

                valueFunctionStationary[i][j][k] = termCost;
                optimalControlsStationary[i][j][k] = 0;
            }
        }
    }

    //--------------------------------Stationary PDE solve-------------------------------------------------

    initControls(gControlGridDim); //init control grid
    stationaryGSSolver(&valueFunctionStationary, &optimalControlsStationary); //SOLVE USING MGS
    rowMajorWriteToFile(&valueFunctionStationary, &optimalControlsStationary);

    gSolvingGreenOnly = false;

    writeToFile(&valueFunctionStationary, &optimalControlsStationary); //write results to file

    //------------------------------Trajectory tracing-------------------------------------------------------

    double tau = gDT;
    gControlInterp = true;
    gStationaryTraj = true;
    gInStationarySolve = true;
    if(gPareto){
        gSelectedTraj = true;
    }

    struct sOptimalTrajectory optTraj = trajectoryTracer(&valueFunctionStationary, &optimalControlsStationary, &valueFunctionStationary, &optimalControlsStationary, &valueFunctionStationary, &optimalControlsStationary, 80, 0, 0, tau);

    cout << "Total discomfort: " << gDiscomfortCost << "\n";
    cout << "Total fuel cost: " << gFuelCost << "\n";
    cout << "Total time: " << gTimeCost << "\n";

    if(gPareto){
        gSelectedTraj = false;
    }

    gControlInterp = false;
    vector<double> optPos = optTraj.positionVals;
    vector<double> optVel = optTraj.velocityVals;
    vector<double> optTime = optTraj.timeVals;
    vector<double> optControls = optTraj.optimalControlVals;
    vector<double> analyticalOCs = optTraj.analyticalOCVals;
    short int loopEndingIdx = optTraj.loopIdx;

    if(!gPareto){
        trajWriteToFile(optPos, optVel, optTime, optControls, loopEndingIdx);
    }

}

/*==============================================================================
 * YR Problem
 *============================================================================*/
void yellowRedProblem(double aTrajStartPos, double aTrajStartVel, double aTrajStartTime){
    /*
     Function: yellowRedProblem

     Arguments: starting (d,v,t) for optimal trajectory tracing.

     Purpose: Solves the red-green OR yellow-red-green problem using the
     two-pass method described in the Optimal Driving manuscript as appropriate.
     Traces optimal trajectory from supplied starting point.

     */

    initGrid();

    short int tPosMaxIdx = (gDNum / 2) + 1;

    //Step 1: Set up d x v matrices to hold solution
    multiarray valueFunctionStationary(boost::extents[gDNum+1][gVNum+1][2]); //array for values of value function at each point in domain
    multiarray optimalControlsStationary(boost::extents[gDNum+1][gVNum+1][2]); //array for acceleration at each point in domain

    multiarray valueFunction(boost::extents[gDNum+1][gVNum+1][gNt + 1]); //array for values of value function at each point in domain
    multiarray optimalControls(boost::extents[gDNum+1][gVNum+1][gNt + 1]); //array for acceleration at each point in domain

    multiarray valueFunctionFirstPass(boost::extents[gDNum+1][gVNum+1][gNt + 1]); //array for values of value function at each point in domain
    multiarray optimalControlsFirstPass(boost::extents[gDNum+1][gVNum+1][gNt + 1]); //array for acceleration at each point in domain

    //Init array for stationary solve
    for(short int i = 0; i < gDNum + 1; i++) {     // distance loop
        for(short int j = 0; j < gVNum + 1; j++) { // velocity loop
            for(short int k = 0; k < 2; k++) {   // time loop

                double d = i * gDPos + gDTarget; //shifted for target at gDTarget

                double distFromTarget = d - gDTarget; //target somewhere to left of d = 0
                double termCost = gMu * pow(distFromTarget, 2); //terminal cost

                valueFunctionStationary[i][j][k] = termCost;
                optimalControlsStationary[i][j][k] = 0;
            }
        }
    }

    //-----------------------------------------PDE solve---------------------------------------------------------

    initControls(gControlGridDim); //init control grid

    gSolvingGreenOnly = true;
    gInStationarySolve = true;
    initControls(gControlGridDim); //init control grid
    stationaryGSSolver(&valueFunctionStationary, &optimalControlsStationary);

    gSolvingGreenOnly = false;

    gInTSolve = false;
    gInStationarySolve = false;
    gControlValues.clear();
    initControls(gControlGridDim); //init control grid

    initGrid();

    //Now, initialize top slices of red / yellow phase arrays
    for(short int i = 0; i < gDNum + 1; i++){     // distance loop
        for(short int j = 0; j < gVNum + 1; j++){ // velocity loop
            for(short int k = gTimeEndIdx; k < gNt + 1; k++){   // time loop
                valueFunction[i][j][k] = valueFunctionStationary[i][j][0];
                optimalControls[i][j][k] = optimalControlsStationary[i][j][0];
            }
        }
    }
    gGreenRedGreenProb = true;
    gFirstPass = true;
    //gValueFunctionValsAtTestPoint.clear(); //debugging array
    timeDepSolver(&valueFunction, &optimalControls, &valueFunctionFirstPass, &optimalControlsFirstPass, &valueFunctionFirstPass, &optimalControlsFirstPass, &valueFunctionStationary, &optimalControlsStationary);

    if(!gPareto){
        rowMajorWriteToFile(&valueFunction, &optimalControls); //for terminal condition in uncertain problem
    }

    if(!gInRGSolve){
        gFirstPass = false;
        valueFunctionYellowPhaseResolve(&valueFunction, &optimalControls, &valueFunctionFirstPass, &optimalControlsFirstPass, &valueFunctionStationary, &optimalControlsStationary);
    }

    if(!gPareto){
        writeToFile(&valueFunction, &optimalControls);
        rowMajorWriteToFile(&valueFunction, &optimalControls);
    }

   //-----------------------------------------------Trajectory tracing------------------------------------------------

    vector<double> fullPosition;
    vector<double> fullVelocity;
    vector<double> fullOC;
    vector<double> fullTime;

    double tau = gDT;
    gControlInterp = true;

    //gPareto = true; //1/7

    if(gPareto){
        gSelectedTraj = true;
    }
    struct sOptimalTrajectory optTraj = trajectoryTracer(&valueFunction, &optimalControls, &valueFunctionFirstPass, &optimalControlsFirstPass, &valueFunctionFirstPass, &optimalControlsFirstPass, aTrajStartPos, aTrajStartVel, aTrajStartTime, tau);

    vector<double> optPos = optTraj.positionVals;
    vector<double> optVel = optTraj.velocityVals;
    vector<double> optTime = optTraj.timeVals;
    vector<double> optControls = optTraj.optimalControlVals;
    vector<double> analyticalOCs = optTraj.analyticalOCVals;
    short int loopEndingIdx = optTraj.loopIdx;

    short int fullSize = loopEndingIdx;

    double startingPos = optPos[loopEndingIdx];
    double startingVel = optVel[loopEndingIdx];

    short int yrPhaseSize = optTime.size();

    if(startingPos > gDTarget){
        gStationaryTraj = true;
        gGreenRedGreenProb = false;
        gInStationarySolve = false;
        gUncertainGreenProb = false;

        struct sOptimalTrajectory optTrajStationary = trajectoryTracer(&valueFunction, &optimalControls, &valueFunction, &optimalControls, &valueFunctionFirstPass, &optimalControlsFirstPass, startingPos, startingVel, 0, tau);


        vector<double> optPosStat = optTrajStationary.positionVals;
        vector<double> optVelStat = optTrajStationary.velocityVals;
        vector<double> optTimeStat = optTrajStationary.timeVals;
        vector<double> optControlsStat = optTrajStationary.optimalControlVals;
        vector<double> analyticalOCsStat = optTrajStationary.analyticalOCVals;

        short int stationaryTrajSize = optTimeStat.size();

        for(short int i = 0; i < stationaryTrajSize; i++) {
            double pos = optPosStat[i];
            double vel = optVelStat[i];
            double time = optTimeStat[i] + gTerminalT;
            double oc = optControlsStat[i];

            optControls.push_back(oc);

            if(i > 0){
                optPos.push_back(pos);
                optVel.push_back(vel);
                optTime.push_back(time);
            }
        }

        fullSize = optTime.size();

    }

    cout << "Total cost: " << gTotalCost << "\n";

    if(gPareto){
        gSelectedTraj = false;
    }
    gControlInterp = false;

    if(!gPareto){
        gInStationarySolve = false;
        gGreenRedGreenProb = true;
        trajWriteToFile(optPos, optVel, optTime, optControls, fullSize);
    }
}

/*==============================================================================
 * GYRG Problem
 *============================================================================*/
void uncertainGreenProblem(double aTrajStartPos, double aTrajStartVel, double aTrajStartTime){
    /*
     Function: uncertainGreenProblem

     Arguments: starting (d,v,t) for optimal trajectory tracing.

     Purpose: Solves the HJB equation(s) for the uncertain turning yellow time problem.
     Reads in results from YRG solve to use as the terminal condition for this problem.
     Traces optimal trajectory from supplied starting point.

     Note: MUST have row-major results from YRG solve with the same priority objectives
     (c1,c2,c3) values saved in .txt file in the same folder to read in the terminal
     condition.
     */

    //Step 1: Solve deterministic yellow-red problem
    initGrid();

    short int tPosMaxIdx = (gDNum / 2) + 1;
    //Step 1: Set up d x v matrices to hold solution
    multiarray valueFunctionStationary(boost::extents[gDNum+1][gVNum+1][2]); //array for values of value function at each point in domain
    multiarray optimalControlsStationary(boost::extents[gDNum+1][gVNum+1][2]); //array for acceleration at each point in domain

    multiarray valueFunction(boost::extents[gDNum+1][gVNum+1][gNt + 1]); //array for values of value function at each point in domain
    multiarray optimalControls(boost::extents[gDNum+1][gVNum+1][gNt + 1]); //array for acceleration at each point in domain

    multiarray valueFunctionFirstPass(boost::extents[gDNum+1][gVNum+1][gNt + 1]); //array for values of value function at each point in domain
    multiarray optimalControlsFirstPass(boost::extents[gDNum+1][gVNum+1][gNt + 1]); //array for acceleration at each point in domain

    gFirstPass = true;
    readFromFile(&valueFunctionFirstPass, &optimalControlsFirstPass);
    gFirstPass = false;
    readFromFile(&valueFunction, &optimalControls);
    initControls(gControlGridDim);

    //Uncertain green problem
    gGreenRedGreenProb = false;
    gInUncertainGreenPhase = true;
    gUncertainGreenProb = true;

    gModel1 = true;
    gModel2 = false;

    gPedestrianModel = false;

    initGrid();

    multiarray valueFunctionUncertGreen(boost::extents[gDNum+1][gVNum+1][gNt + 1]); //array for values of value function at each point in domain
    multiarray optimalControlsUncertGreen(boost::extents[gDNum+1][gVNum+1][gNt + 1]); //array for acceleration at each point in domain

    //Now, initialize top slices of red / yellow phase arrays
    for(short int i = 0; i < gDNum + 1; i++){     // distance loop
        for(short int j = 0; j < gVNum + 1; j++){ // velocity loop
            for(short int k = gNt; k >= 0; k--){ // velocity loop
                valueFunctionUncertGreen[i][j][k] = valueFunction[i][j][0];
                optimalControlsUncertGreen[i][j][k] = optimalControls[i][j][0];
            }
        }
    }

    timeDepSolver(&valueFunctionUncertGreen, &optimalControlsUncertGreen, &valueFunction, &optimalControls, &valueFunctionFirstPass, &optimalControlsFirstPass, &valueFunctionStationary, &optimalControlsStationary);

    if(!gPareto){
        writeToFile(&valueFunctionUncertGreen, &optimalControlsUncertGreen);
    }

    double tau = gDT;
    gControlInterp = true;
    gModel1Traj = (gModel1) ? true : false;
    struct sOptimalTrajectory optTraj;

    if(gPareto){
        gSelectedTraj = true;
        optTraj = trajectoryTracer(&valueFunctionUncertGreen, &optimalControlsUncertGreen, &valueFunction, &optimalControls, &valueFunctionFirstPass, &optimalControlsFirstPass, aTrajStartPos, aTrajStartVel, aTrajStartTime, tau);
        gSelectedTraj = false;
    }

    else{
        if(gNumStops == 2){
            model1TwoLengthsTrajectories(&valueFunctionUncertGreen, &optimalControlsUncertGreen, &valueFunction, &optimalControls, &valueFunctionFirstPass, &optimalControlsFirstPass, &valueFunctionStationary, &optimalControlsStationary, aTrajStartPos, aTrajStartVel, aTrajStartTime, tau);
        }

        else if(gNumStops == 3){
            model1ThreeLengthsTrajectories(&valueFunctionUncertGreen, &optimalControlsUncertGreen, &valueFunction, &optimalControls, &valueFunctionFirstPass, &optimalControlsFirstPass, &valueFunctionStationary, &optimalControlsStationary, aTrajStartPos, aTrajStartVel, aTrajStartTime, tau);
        }

        else if(gNumStops == 4){
            model1FourLengthsTrajectories(&valueFunctionUncertGreen, &optimalControlsUncertGreen, &valueFunction, &optimalControls, &valueFunctionFirstPass, &optimalControlsFirstPass, &valueFunctionStationary, &optimalControlsStationary, aTrajStartPos, aTrajStartVel, aTrajStartTime, tau);
        }
    }

    gInUncertainGreenPhase = false;
    gControlInterp = false; //need to reset for each pareto solve
     //Full traj
     gModel1 = false;
     gModel2 = false;

}

/*==============================================================================
 * Run Optimal Driving Examples
 *============================================================================*/
void runOptimalDrivingExamples(int aProbSelectionIdx){
  /*
   Function: runOptimalDrivingExamples

   Arguments: integer corresponding to example number
   See Optimal Driving manuscript for the example numbers.
   aProbSelectionIdx is passed to the inerface as a command-line argument.

   Purpose: Solves PDE(s) and obtains all results files needed to generate
   the plots associated with the selected example.
   */

  //Example 1: Green only (stationary)
  if(aProbSelectionIdx == 1){
      //Run the trajectories first
      gC1 = 0.025;
      gC2 = 0.025;
      gC3 = 0.95;

      gInExample1TrajectoryComputation = true;
      greenLightOnlyProblem(80, 0, 0);
      gInExample1TrajectoryComputation = false;
      setFrameworkSettings();

      gC1 = 1.0/3.0;
      gC2 = 1.0/3.0;
      gC3 = 1.0/3.0;

      gPareto = true;
      gStationaryPareto = true;
      paretoGrid(40, 40, gC1, gC2, gC3, 80, 0, 0);

      gStationaryPareto = true;
      gJ3ConstrainedCurve1 = true;
      j3ConstrainedParetoFront(40, 40, 45, 80, 0, 0);
      gJ3ConstrainedCurve1 = false;

      gJ3ConstrainedCurve2 = true;
      j3ConstrainedParetoFront(40, 40, 35, 80, 0, 0);
      gJ3ConstrainedCurve2 = false;

      gJ3ConstrainedCurve3 = true;
      j3ConstrainedParetoFront(40, 40, 25, 80, 0, 0);
      gJ3ConstrainedCurve3 = false;
  }

  //Example 2: Red-Green problem
  else if (aProbSelectionIdx == 2){
      gInRGSolve = true;
      cout << "Computing first trajectory..." << "\n";
      yellowRedProblem(80, 15, 0);
      setFrameworkSettings();
      gInRGSolve = true;
      cout << "Computing second trajectory..." << "\n";
      yellowRedProblem(80, 0, 0);

  }

  //Example 3: Yellow-Red-Green problem
  else if (aProbSelectionIdx == 3){
      yellowRedProblem(43, 10, 0);
      cout << "Computing first trajectory..." << "\n";
      cout << "Done." << "\n";
      setFrameworkSettings();
      cout << "Computing second trajectory..." << "\n";
      yellowRedProblem(48, 10, 0);
      cout << "Done.";
  }

  //Example 4: Uncertain TY, 2 possible turning yellow times
  else if (aProbSelectionIdx == 4){
      gNumStops = 2;
      //c1 = c2 = c3 solve first
      cout << "Computing first example..." << "\n";
      yellowRedProblem(94, 0.85, 0);
      setFrameworkSettings();
      gInEqualCEqualP2lExample = true;
      uncertainGreenProblem(94, 0.85, 0);

      cout << "Computing second example..." << "\n";
      setFrameworkSettings();
      gInEqualCEqualP2lExample = false;
      uncertainGreenProblem(94, 0.85, 0);

      //c2 prioritized computation
      cout << "Computing third example..." << "\n";
      setFrameworkSettings();
      gNumStops = 2;
      gInC2Priority2lExample = true;
      yellowRedProblem(94, 0.85, 0);
      setFrameworkSettings();
      gInC2Priority2lExample = true;
      uncertainGreenProblem(94, 0.85, 0);

  }

  //Example 5: Uncertain TY, 3 possible turning yellow times
  else if (aProbSelectionIdx == 5){
      gNumStops = 3;
      yellowRedProblem(68, 5, 0);
      setFrameworkSettings();
      gInEqualCEqualP2lExample = true;
      uncertainGreenProblem(68, 5, 0);
  }

}
