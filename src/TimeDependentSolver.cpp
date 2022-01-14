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
 * File: TimeDependentSolver.cpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains the function responsible for executing the
 * time-dependent HJB solve during the red phase, yellow phase, and the
 * initial uncertain green phase. This solver solves the time-dependent HJB
 * equations presented in the manuscript backwards in time using a
 * Semi-Lagrangian discretization.
 *
 * This file also contains the function responsible for executing the second-pass
 * sweep during the yellow phase in the YRG problem to account for any discontinuity
 * present in the problem domain resulting from the driver's objective preferences
 * and the tradeoff between having the ability to speed through the intersection
 * before the light turns red or come to a complete stop before the red light.
 *============================================================================*/

#include "TimeDependentSolver.hpp"

//---------------Libraries----------------------------------------------------
#include <iostream>
#include <fstream>
#include <cmath>

//------------------Project-specific header files----------------------------
#include "boost/multi_array.hpp"
#include "variables.h"
#include "ProblemFunctions.hpp"
#include "InitFunctions.hpp"
#include "StationarySolver.hpp"
#include "CoordinateFunctions.hpp"
#include "ProbabilityFunctions.hpp"
#include "HelperFunctions.hpp"

/*==============================================================================
* Time Dependent HJB Solver
*============================================================================*/
void timeDepSolver(multiarray *aU, multiarray *aOC, multiarray *aUYRG, multiarray *aOCYRG, multiarray *aUYRGFirstPass, multiarray *aOCYRGFirstPass, multiarray *aStationarySoln, multiarray *aStationaryOC){
    /*
     Function: timeDepSolver

     Arguments: arrays to hold the value function and optimal controls

     Purpose: Solves the time dependent HJB equation(s) associated with the selected
     example problem backwards in time using a Semi-Lagrangian discretization. The solver
     saves a value function value and an optimal control value at each node in the supplied
     grid. This function relies heavily on the controlLoop, isAllowed, prescribedValueFunction, conditionalProb, isOnParabPiece, and isOnMaxAccelBD functions.
     */
    short int turningRedIdx = short(int(round(gTR / gDT))); //time index of turning red time
    short int idxAtLight = (gDNum / 2) + 1 ; //position index at d_l
    assert(turningRedIdx * gDT == gTR);
    for(short int k = gNt-1; k >= gTimeEndIdx; k--) {//time loop
        double currentTime = k * gDT; //current physical time
        gTimeIndex = k;

        if(!gPareto){
            cout << "k: " << k << "\n"; //countdown ticker to track solving progress
        }

        double pStar = conditionalProb(currentTime + gDT); //probability of light change at t+\Delta t
        pStar = (gInUncertainGreenPhase) ? pStar : 0; //set pStar = 0 if in deterministic phase

        gProbVals.push_back(pStar); //debugging - array of probability values

        if ((currentTime < gTG) && (gGreenRedGreenProb)){
            gHasBoundaryOn = true; //turn on pw boundaries once we get to red phase interval
        }

        else if (gInUncertainGreenPhase){
            gHasBoundaryOn = false; //turn off pw boundaries during uncertain green phase
        }

        for(short int i = idxAtLight; i < gDNum + 1; i++) {// distance loop
            double currentPosition = i * gDPos + gDTarget; //current physical pos
            gPosIndex = i;
            for(short int j = 0; j < gVNum + 1; j++) {// velocity loop
                double currentVelocity = j * gDVel; //current physical vel
                gVelIndex = j;

                double bestVal = gInfty; //place holder for value function val
                double optimalControl = gAccelInfty; //placeholder for optimal control val

                gTau = gDT; //set SL timestep = \Delta t

                bool startsInLegalRegion = isAllowed(currentVelocity, currentPosition, currentTime); //determine if node is in allowable region

                double vfIsPrescribed = prescribedValueFunction(aU, aOC, currentPosition, currentVelocity, currentTime, gDT, 0); //determine if value function known already at node

                double dParab = pow(currentVelocity, 2) / (2 * gBeta); //location of parab of last resort

                //set value function if it is known already
                if(vfIsPrescribed >= 0){
                    bestVal = vfIsPrescribed;
                    if(gHasBoundaryOn){
                        bool isOnParab = isOnParabPiece(currentPosition, currentVelocity, currentTime);
                        bool onMaxAccelBd = isOnMaxAccelBD(currentPosition, currentVelocity, currentTime);
                        if ((isOnParab) || (currentPosition >= dParab)){
                            assert(currentVelocity != 0);
                            optimalControl = -gAlpha; //max braking
                        }

                        else if (onMaxAccelBd){
                            optimalControl = gBeta; //max accel
                        }

                        else{
                            //on target
                            optimalControl = 0;
                        }
                    }

                    else{
                        optimalControl = 0;
                    }
                }
                //point in allowable region at t but need to interpolate to get value function and OC
                else if (startsInLegalRegion){
                    struct sOptimalValues optVals = controlLoop(aU, aOC, aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, currentPosition, currentVelocity, currentTime, k, gDT, pStar); //call control loop to determine value function and OC

                    //set value function and oc
                    bestVal = optVals.valueFn;
                    optimalControl = optVals.oc;
                }

                //save value function and oc to multiarray before continuing
                (*aU)[i][j][k] = bestVal;
                (*aOC)[i][j][k] = optimalControl;

            }//end vel loop
        }//end pos loop
    }//end k loop
    cout << "Through first pass solve" << "\n";
}

/*==============================================================================
 * Value function corrector for [TY, TR)
 *============================================================================*/
void valueFunctionYellowPhaseResolve(multiarray *aValueFunction, multiarray *aOptimalControls, multiarray *aValueFunctionFirstPass, multiarray *aOptimalControlsFirstPass, multiarray *aStationarySoln, multiarray *aStationaryOC){
    /*
     Function: valueFunctionYellowPhaseResolve

     Arguments: arrays to hold the value function and optimal controls. When this function
     is called, the value function and OC arrays should have already been pre-populated
     using the results from the parabolic bd-only first pass.

     Purpose: Correction method applied to determine the value function during the yellow phase. This method passes backwards in time from the start of the red phase, down to the start of the yellow phase, and resets the value function and optimal control value at (i,j,k) if the value computed after taking any discontinuity into account is smaller than the first pass value. Only values to the RIGHT of the light at d = 0 and to the LEFT or ON the potential discontinuity line are candidates for resetting. This function relies heavily on the controlLoop, isAllowed, prescribedValueFunction, conditionalProb, isOnParabPiece, and isOnMaxAccelBD functions.
     */

    short int turningRedIdx = round(gTR / gDT); //time index of red light
    short int idxAtLight = gDNum / 2; //position index at traffic light at d = 0
    for(short int k = turningRedIdx - 1; k >= 0; k--) {//time loop
        double currentTime = k * gDT; //physical time
        gTimeIndex = k;

        if(!gPareto){
            //cout << "k: " << k << "\n"; //ticker to track progress
        }

        double pStar = conditionalProb(currentTime + gDT);
        pStar = (gInUncertainGreenPhase) ? pStar : 0;

        if ((currentTime < gTG) && (gGreenRedGreenProb)){
            gHasBoundaryOn = true; //turn on pw boundaries once we get to red phase interval
        }

        for(short int i = idxAtLight; i < gDNum + 1; i++) {// distance loop
            double currentPosition = i * gDPos + gDTarget;
            gPosIndex = i;
            for(short int j = 0; j < gVNum + 1; j++) {// velocity loop
                //Compare and set value function
                double currentVelocity = j * gDVel;
                gVelIndex = j;

                double valueFunctionValueWithFullParabBd = (*aValueFunction)[i][j][k];

                double valueFunctionValue = gInfty;
                double optimalControlValue = gAccelInfty;

                gTau = gDT; //tau reset
                gFirstPass = false;

                //---------------------Cases considered based on where we are in the process-------------
                bool startsInLegalRegion = isAllowed(currentVelocity, currentPosition, currentTime);

                double vfIsPrescribed = prescribedValueFunction(aValueFunction, aOptimalControls, currentPosition, currentVelocity, currentTime, gDT, 0);
                double dParab = pow(currentVelocity, 2) / (2 * gBeta);
                double maxAccelLinePos = maxAccelLine(currentVelocity, currentTime);

                if(vfIsPrescribed >= 0){
                    valueFunctionValue = vfIsPrescribed;
                    if(gHasBoundaryOn){
                        bool isOnParab = isOnParabPiece(currentPosition, currentVelocity, currentTime);
                        bool onMaxAccelBd = isOnMaxAccelBD(currentPosition, currentVelocity, currentTime);
                        if ((isOnParab) || (currentPosition >= dParab)){
                            assert(currentVelocity != 0);
                            optimalControlValue = -gAlpha; //max braking
                        }

                        else if (onMaxAccelBd){
                            optimalControlValue = gBeta;
                        }

                        else{
                            //on target
                            optimalControlValue = 0;
                        }
                    }

                    else{
                        optimalControlValue = 0;
                    }
                }

                else if (startsInLegalRegion){
                    //cout << "legal start" << "\n";
                    struct sOptimalValues optVals = controlLoop(aValueFunction, aOptimalControls, aValueFunctionFirstPass, aOptimalControlsFirstPass, aValueFunctionFirstPass, aOptimalControlsFirstPass, currentPosition, currentVelocity, currentTime, k, gDT, pStar);

                    valueFunctionValue = optVals.valueFn;
                    optimalControlValue = optVals.oc;

                    /*if (valueFunctionValue >= gInfty){
                        cout << "i: " << gPosIndex << " j: " << gVelIndex << " k: " << gTimeIndex << "\n";
                    }*/ //debugging

                    assert(valueFunctionValue < gInfty);
                }

                //Value comparison: Reset if smaller than first pass value
                if((valueFunctionValue < valueFunctionValueWithFullParabBd) && (currentPosition < maxAccelLinePos)){
                    (*aValueFunction)[i][j][k] = valueFunctionValue;
                    (*aOptimalControls)[i][j][k] = optimalControlValue;
                }
            }
        }
    }

    cout << "Through second pass / correction solve" << "\n";
}
