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
 * File: InitFunctions.cpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains the functions utilized in initializing
 * the grids, parameters, and variables used throughout Optimal Driving.
 *============================================================================*/

//---------------------Libraries----------------------------------------------
#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>
#include <utility>

//---------------------Project-specific header files----------------------------
#include "variables.h"
#include "InitFunctions.hpp"

using namespace std;

/*==============================================================================
 * Set framework settings
 *============================================================================*/
void setFrameworkSettings(){
    /*
     Function: setFrameworkSettings()

     Purpose: Sets all of the boolean switches used throughout the Optimal Driving Framework.

     See main.cpp for complete descriptions of their purposes.

     NOTE: These switches should NOT be manually toggled once set.
     */
    //=========================================GRID REF CONTROL==================================================
    gGridRef = true;
    gGridRefFactor = 4;

    //========================================PROBLEM SELECTION SWITCHES==========================================
    gInStationarySolve = false;
    gSolvingGreenOnly = false;
    gInGS = false;
    gModifiedGS = true;
    gInGJ = false;

    gGreenRedGreenProb = false; //on for grg problem
    gInTSolve = false;
    gInRGSolve = false;

    gUncertainGreenProb = false;
    gInUncertainGreenPhase = false; //boolean to indicate if we are solving the PDE in the uncertain initial green phase
    gInTerminalPenalty = false;

    //=========================================UNCERTAINTY SWITCHES / INIT========================================
    gModel1 = false;
    gModel2 = false;
    gPedestrianModel = false;
    gInEqualCEqualP2lExample = false;
    gInC2Priority2lExample = false;

    gLightTimes.clear(); //reset turning yellow times array
    gLightProbs.clear(); //reset turning yellow probabilities array
    //========================================INTERPOLATION SWITCHES==============================================
    gHasBoundaryOn = false;
    gFirstPass = false; //indicates if in first sweep of YRG problem
    gInMaxAccelInterp = false; //indicates if in interpolation function along max accel line
    gInParabInterp = false;//indicates if in parabolic interpolation function

    //=========================================TRAJECTORY TRACER SWITCHES=========================================
    gInTracer = false; //indicates if in trajectory tracer
    gControlInterp = false; //indicates if using control interpolation in tracer
    gStationaryTraj = false; //indicates if tracing trajectory for the stationary problem

    gModel1Traj = false; //indicates if model 1 uncertainty used and model 1 trajectory is being traced
    gModel1TrajDetPortion = false; //indicates if tracing deterministic portion of a model 1 trajectory
    gFirstTraj = false;
    gSecondTraj = false;
    gThirdTraj = false;
    gFourthTraj = false;

    //======================================PARETO SWITCHES=======================================================
    gPareto = false; //indicates if we are finding Pareto front
    gSelectedTraj = false; //indicates trajectory tracing in Pareto framework
    gParetoCostComputation = false; //indicates if we are computing cost along trajectory in Pareto
    gStationaryPareto = false;
    gYRGPareto = false;
    gUncertainPareto = false;
    gInConstrainedPareto = false;
    gJ3ConstrainedCurve1 = false;
    gJ3ConstrainedCurve2 = false;
    gJ3ConstrainedCurve3 = false;

    gFuelCost = 0;
    gDiscomfortCost = 0;
    gTimeCost = 0;
    gTotalCost = 0;
}

/*==============================================================================
 * Init grid
 *============================================================================*/
void initGrid(){
        /*
         Function: initGrid()

         Purpose: Initializes computational grid with selected paramater type
         according to which problem we are in.

         Artificial parameters (not used) initializes a grid where the max position
         and velocity are set to one, min/max allowed accel is set to -1 and 1
         respectively, and the target sits at position -1.

         Realistic parameters (used throughout all examples in the manuscript)
         initalizes a grid where the bounds and light durations are set to reflect
         a real-world driving scenario.
         */
    //---------------Light lengths and yellow light info------------------------
    double greenDuration;

    if (gArtificial){
        gTY1 = 0.2; //first turning yellow time
        gTYMax = 0.6; //second turning yellow time
        gTYP1 = 0.5; //probability of TY1
        gTYP2 = 0.5; //probability of TY2
        greenDuration = 4; //how long the green light length post-red phase should be.
        gYellowDuration = 0.5; //length of yellow light (seconds)
        gRedDuration = 0.6; //length of red light (seconds)

        gVMax = 1; //max velocity (meters / second)
        gDMax = 1; //max distance (meters)
        gDTarget = -1; //max distance (meters)
        gAlpha = 1; //max decel (meters / second / second)
        gBeta = 1; //max accel
        gAccelInfty = -gAlpha - 3;

        //--------------------------d,v,t grid-------------------------------------

        gN = 100;
        gVNum = gN; //number of velocity intervals to take

        cout << "Nv: " << gVNum << "\n";
        gDVel = gVMax / gN; //delta v

        gAccelInfty = -gAlpha - 1;
    }//end artificial

    //-----------------------Init Realistic params-----------------------------------
    else{
        //Reset c values if in Example 4 when c2 is prioritized over the others
        if(gInC2Priority2lExample){
            gC1 = 0.15;
            gC2 = 0.75;
            gC3 = 0.1;

            gTYP1 = 0.5; //probability of TY1
            gTYP2 = 0.5; //probability of TY2
        }

        //Reset yellow light change probabilities if in Example 4 when T_1 is most likely
        if((!gInEqualCEqualP2lExample) && (!gInC2Priority2lExample) && (gNumStops == 2)){
            gTYP1 = 0.95;
            gTYP2 = 0.05;
        }

        //Reset yellow light change probabilities in Example 4
        else if (gNumStops == 2){
            gTYP1 = 0.5; //probability of TY1
            gTYP2 = 0.5; //probability of TY2
        }

        //-----------------------------Parameter Values---------------------------------------
        gTY1 = 2; //first turning yellow time
        gTYMax = 6; //second turning yellow time //To do: needs to be longer! 30 seconds too short for green light
        gYellowDuration = 3; //length of yellow light (seconds)
        gYellowDuration = (gDebugControlChatter) ? gTYMax - gTY1 + gYellowDuration : gYellowDuration; //extend horizon for the debuging

        cout << "gYellowDuration:" << gYellowDuration << "\n";

        gRedDuration = 60; //length of red light (seconds)
        greenDuration = 120; //length of final green light (not in use)

        gVMax = 20.11; //max velocity (meters / second)
        gDMax = 100; //max distance (meters)
        gDTarget = -100; //max distance (meters)
        gAlpha = 3.8; //max decel (meters / second / second)
        gBeta = 3.8; //max accel
        gAccelInfty = -gAlpha - 3;

        //-------------------------------Velocity grid--------------------------------------------------
        gN = 45; //Number of velocity INTERVALS

    }//end realistic else

    //======================================Grid Creation=============================================================

    //Example 4, two possible turning yellow times
    if((gNumStops == 2) && (gUncertainGreenProb)){

        gTY1 = 2; //first turning yellow time
        gTY2 = 6; //second turning yellow time
        gTYMax = gTY2; //second turning yellow time

        //Set probs if in example where T1 prioritized
        if((!gInEqualCEqualP2lExample) && (!gInC2Priority2lExample)){
            gTYP1 = 0.95;
            gTYP2 = 0.05;
        }

        //Set probs if in equal prob examples
        else{
            gTYP1 = 0.5; //probability of TY1
            gTYP2 = 0.5; //probability of TY2
        }

        //Set objective priorities if in example where c2 prioritized
        if(gInC2Priority2lExample){
            gC1 = 0.15;
            gC2 = 0.75;
            gC3 = 0.1;
        }

        else{
            gC1 = 1.0/3.0;
            gC2 = 1.0/3.0;
            gC3 = 1.0/3.0;
        }

        //Save possible turning yellow times and associated probabilities
        gLightTimes.push_back(gTY1);
        gLightTimes.push_back(gTYMax);

        gLightProbs.push_back(gTYP1);
        gLightProbs.push_back(gTYP2);

    }

    //Example 5, three possible turning yellow times
    else if((gNumStops == 3) && (gUncertainGreenProb)){
        gTY1 = 2; //first turning yellow time
        gTY2 = 4; //second turning yellow time
        gTY3 = 6; //third turning yellow time
        gTYMax = gTY3; //max turning yellow time - longest possible initial uncertain green phase duration
        gTYP1 = 0.25; //probability of TY1
        gTYP2 = 0.25; //probability of TY2
        gTYP3 = 0.5; //probability of TY3

        //Save light times and probabilities
        gLightTimes.push_back(gTY1);
        gLightTimes.push_back(gTY2);
        gLightTimes.push_back(gTY3);

        gLightProbs.push_back(gTYP1);
        gLightProbs.push_back(gTYP2);
        gLightProbs.push_back(gTYP3);
    }

    //Four possible turning yellow times
    else if ((gNumStops == 4) && (gUncertainGreenProb)){
        //Light Lengths
        gTY1 = 2;
        gTY2 = 4;
        gTY3 = 5;
        gTY4 = 6;

        gTYMax = gTY4;

        //Probs:
        gTYP1 = 0.25;
        gTYP2 = 0.25;
        gTYP3 = 0.25;
        gTYP4 = 0.25;

        //Saving to array
        gLightTimes.push_back(gTY1);
        gLightTimes.push_back(gTY2);
        gLightTimes.push_back(gTY3);
        gLightTimes.push_back(gTY4);

        gLightProbs.push_back(gTYP1);
        gLightProbs.push_back(gTYP2);
        gLightProbs.push_back(gTYP3);
        gLightProbs.push_back(gTYP4);
    }

    gTYMax = (gInUncertainGreenPhase) ? gTYMax : 0; //set the start of the deterministic yellow phase

    gTR = gTYMax + gYellowDuration; //set turning red time based on start of yellow
    gTG = gTR + gRedDuration; //set turning green time based on turning red time and duration of red

    gTerminalT = gYellowDuration + gRedDuration; //set terminal time of the deterministic, time dependent process

    //If solving red-green only, make t = 0 the turning red time
    if(gInRGSolve){
        gTR = 0;
        gTG = gTR + gRedDuration;
        gTerminalT = gTG;
    }

    gTerminalT = (gInUncertainGreenPhase) ? gTYMax : gTerminalT; //set terminal time based on if we are in uncertain stage or not

    //------------------------------Velocity grid--------------------------------------------

    gVNum = gN; //number of velocity intervals
    gDVel = gVMax / gVNum; //Delta v

    cout << "Nv: " << gVNum << "\n";
    cout << "Delta v: " << gDVel << "\n";

    //-------------------------------Position Grid---------------------------------------
    double maxAbsoluteAccel = max(gAlpha, gBeta); //maximum of absolute value of min/max allowed accel
    double provisionalDeltaT = gDVel / maxAbsoluteAccel; //set provisional \Delta t based on max absolute allowed accel

    double provisionalDeltaD = (gVMax * provisionalDeltaT); //set provisional \Delta d based on largest amount allowed to traverse in \Delta t sec
    gDNum = ceil(fabs(gDMax - gDTarget) / provisionalDeltaD); //number of position intervals
    gDNum = (gDNum % 2 == 0) ? gDNum : gDNum + 1; //ensure number of pos intervals even so that d = 0 a gridline
    gDPos = fabs(gDMax - gDTarget) / gDNum; //Delta d

    cout << "Nd: " << gDNum << "\n";
    cout << "Delta d: " << gDPos << "\n";
    cout << "0 idx conversion: " << (gDNum / 2) * gDPos + gDTarget << "\n";
    assert((gDNum / 2) * gDPos + gDTarget == 0);

    //----------------------------------------------Time Grid-----------------------------------------------
    gDT = min(provisionalDeltaT, gDPos / gVMax);
    gTau = gDT; //set tau to be \Delta t to start
    gNt = ceil(gTerminalT / gDT); //number of time intervals

    gDT = 0.1; //setting gDT = 0.1 to keep each second on a time gridline
    gTau = gDT;
    gNt = ceil(gTerminalT / gDT); //set number of time intervals whether we are planning in deterministic or uncertain phase

    cout << "Nt: " << gNt << "\n";
    cout << "Delta t: " << gDT << "\n";

    //set stationary time intervals and terminal time (Note: not relevant in PDE solve, only set high in trajectory tracer to avoid issues with not tracing the entire trajectory due to timestep restrictions put in place for tracing during time dep parts of the problem)
    if((gInStationarySolve)){
        gNt = 1000000;
        gTerminalT = 10000000;

        if(gInExample1TrajectoryComputation){
            gC1 = 0.025;
            gC2 = 0.025;
            gC3 = 0.95;
        }
    }

    //-----------------------------------------Grid refinement---------------------------------------
    if (gGridRef){
        gridRef(gDPos, gDVel, gDT, gGridRefFactor);
    }
}

/*==============================================================================
 * Control Array Init
 *============================================================================*/
void initControls(int aControlArrayDim){
    /*
     Function: initControls()

     Purpose: Initializes grid of aControlArrayDim control values to use in initial
     grid search for the optimal control in the PDE solves before using Golden Section
     Search (GSS) to refine the optimal control approximation.
     */

    gControlValues.clear(); //reset control value array

    //Go through -alpha, 0 to set equally spaced control values
    double accelMin = -gAlpha;
    double accelMax = gBeta;

    //Set the ends of the interval
    gControlValues.push_back(accelMin);
    double controlSpacing = fabs(accelMax - accelMin) / (gNumControls - 1);
    cout << controlSpacing << "\n";
    gControlSpacing = controlSpacing;

    //Loop over controls to set interior values
    for(short int i = 1; i < gNumControls; i++) {
        double controlVal = accelMin + (i * controlSpacing);
        controlVal = (fabs(controlVal) <= 1e-14) ? 0 : controlVal;
        gControlValues.push_back(controlVal);
    }//end loop over controls

    //Modification - Only consider nonnegative control values in stationary problem
    if(gInStationarySolve){
        for(short int i = 0; i < gNumControls; i++) {
            double controlVal = gControlValues[i];
            if(controlVal < 0){
                gControlValues[i] = 0;
            }
        }//end loop over controls
    }

    //Print control values to check to make sure values are correct, if desired
    /*if(!gPareto){
        for(short int i = 0; i < gNumControls; i++) {
            double controlVal = gControlValues[i];
            cout << controlVal << "\n";
        }//end loop over controls
    }*/
}

/*==============================================================================
 * Array initialization
 *============================================================================*/
//void initArray(multiarray *aValueFnArray, multiarray *aOptimalControlArray, const double aOptControl){
    /*
     * Function: initArray()
     *
     * Purpose: Populates a boost multiarray for the value function and optimal controls
     * with the appropriate
     */

    /*for(short int i = 0; i < gDNum + 1; i++) {     // distance loop
        for(short int j = 0; j < gVNum + 1; j++) { // velocity loop
            for(short int k = 0; k < gNt+1; k++) {   // time loop

                //Distance, velocity, and time values
                double d = i * gDPos + gDTarget; //shifted for target at -1
                double v = j * gDVel;
                double t = k * gDT;

                //Setting the values at the final timestep
                if (k==gNt){
                    //TODO: Fix for stationary
                    //(*aValueFnArray)[i][j][k] = terminalCost(aValueFnArray, d, v, t, gMu, gDTarget);
                    //(*aOptimalControlArray)[i][j][k] = aOptControl;
                }

                //Setting the values at all other timesteps (this bit doesn't matter since it all gets rewritten in DP loop)
                else{
                    (*aValueFnArray)[i][j][k] = 0;
                    (*aOptimalControlArray)[i][j][k] = aOptControl;
                }
            }
        }
    }
}*/


/*==============================================================================
 * Grid refinement
 *============================================================================*/
void gridRef (double aDeltaD, double aDeltaV, double aDeltaT, double aRefinementFactor){
    /*
     Function: gridRef()

     Purpose: Refines original grid by aRefinementFactor and updates original aDeltaD,
     aDeltaV, and aDeltaT accordingly.
     */

    //--------------------------Refine \Delta d, \Delta v, \Delta t-------------------------------------
    double refinementFrac = 1 / aRefinementFactor;
    cout << refinementFrac << "\n";
    gDVel = aDeltaV * refinementFrac;
    gDPos = aDeltaD * refinementFrac;
    gDT = aDeltaT * refinementFrac;
    gTau = gDT;
    //--------------------------------Adjust v,d,t interval numbers------------------------------------
    gVNum = gVMax / gDVel;
    gN = gVNum;
    gDNum = fabs(gDMax - gDTarget) / gDPos;
    gNt = ceil(gTerminalT / gDT);

    //----------------------------------Report updated params-----------------------------------------
    cout << "Nv: " << gVNum << "\n";
    cout << "Delta v: " << gDVel << "\n";

    cout << "Nd: " << gDNum << "\n";
    cout << "Delta d: " << gDPos << "\n";
    cout << "0 idx conversion: " << (gDNum / 2) * gDPos + gDTarget << "\n";

    cout << "Nt: " << gNt << "\n";
    cout << "Delta t: " << gDT << "\n";

    assert((gDNum / 2) * gDPos + gDTarget == 0); //ensure that d = 0 still a position gridline

}
