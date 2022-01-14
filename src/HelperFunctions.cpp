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
 * File: HelperFunctions.cpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains all helper functions used in Optimal Driving,
 which perform tasks such as optimal control grid search loop, Golden Section Search,
 running cost integration, boundary costs, trajectory tracing, writing output to file,
 reading data inputs from files, etc.
 *============================================================================*/

#include "HelperFunctions.hpp"

//-------------------Libraries----------------------------------------------
#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <cmath>

//------------------Project specific header files--------------------------
#include "variables.h"
#include "InitFunctions.hpp"
#include "CoordinateFunctions.hpp"
#include "InterpolationFunctions.hpp"
#include "ProblemFunctions.hpp"
#include "ProbabilityFunctions.hpp"
#include "boost/multi_array.hpp"

using namespace std;

typedef boost::multi_array<double, 3> multiarray;

/*==============================================================================
 * Control Grid Search and Optimal Control Determination
 *============================================================================*/
struct sOptimalValues controlLoop(multiarray *aU, multiarray *aOC, multiarray *aUYRG, multiarray *aOCYRG, multiarray *aUFirstPassYRG, multiarray *aOCFirstPassYRG, double aCurrentPos, double aCurrentVel, double aCurrentTime, short int aCurrentTimeIndex, double aTau, double aProbOfLightChange){
    /*
     Function: controlLoop()

     Arguments: value function and OC arrays, current physical position, velocity, and time, current timeslice index,
     SL timestep (tau), probability of light change (only relevant during uncertain phase)

     Purpose: Runs loop over control values to determine optimal control and value function at (i,j,k).
     After the optimal control grid search is complete, it calls Golden Section Search to further refine
     the optimal control approximation at (i,j,k).

     Returns: Struct containing value function value and optimal control value at node (i,j,k). Also
     contains the dNext, vNext, and tNext values computed using the optimal control value determined at that point.
     These dNext, vNext, tNext values are relevant when trajectory tracing using the PDE solution.
     */

    sOptimalValues S;

    short int currentTimeIdx = aCurrentTimeIndex;

    double dNextOpt;
    double vNextOpt;
    double tNextOpt;

    double bestValueFnVal = gInfty; //for loop comparison
    double bestOCVal = gAccelInfty; //for loop comparison

    for(int controlIndex = 0; controlIndex < gNumControls; controlIndex++){

        double currentValue; //value function value placeholder
        double currentControl = gControlValues[controlIndex]; //considered control value

        struct sNextCoordinates nextCoords = tauRefinement(aCurrentPos, aCurrentVel, aCurrentTime, currentControl, aTau);
        double dNext = nextCoords.nextPos;
        double vNext = nextCoords.nextVel;
        double tNext = nextCoords.nextTime;
        gTau = nextCoords.tauValue;

        short int nextTimeIdx = currentTimeIdx + 1;

        if (gInTracer){
            nextTimeIdx = round(tNext / gDT);
        }

        //turn bd on during YR phases
        if ((gGreenRedGreenProb)){
            gHasBoundaryOn = (tNext >= gTG) ? false : true;
        }

        //check to make sure pw boundary off during uncertain green phase
        if(gInUncertainGreenPhase){
            assert(!gHasBoundaryOn);
        }

        bool isOKAtTNext = isAllowed(vNext, dNext, tNext); //determine if next point in allowable region

        double nextVal;
        double nextValParab = 0;

        double delta = terminalCost(aU, aOC, aUFirstPassYRG, aOCFirstPassYRG, dNext, vNext, tNext, aCurrentPos, aCurrentVel, aCurrentTime, currentControl, gMu, gDTarget); //terminal cost - only relevant during uncertain phase

        double vfIsPrescribed = prescribedValueFunction(aU, aOC, dNext, vNext, tNext, aTau, currentControl); //determine if point on a boundary, set value function accordingly

        //Make sure that no boundaries and such are on after light turns green
        if (tNext >= gTG){
            assert(vfIsPrescribed == -10);
        }
        //trajectory tracing with smaller tau - not used in any examples in the manuscript
        if ((gInTracer) && (aTau < gDT)){
            nextVal = dvtInterpolation(aU, aOC, aUYRG, aOCYRG, vNext, dNext, aCurrentPos, aCurrentVel, aCurrentTime, tNext, currentControl);

            double runningCostInt = runningCostIntegral(currentControl, aCurrentTime, aTau, gGamma, gC1, gC2, gC3);

            currentValue = runningCostInt + (1 - aProbOfLightChange) * nextVal + aProbOfLightChange * delta;
        }

        //determine value function and OC in PDE solve
        else{
            if (!isOKAtTNext){
                nextVal = gInfty;
                nextValParab = gInfty;
            }

            else if (vfIsPrescribed >= 0){
                nextVal = vfIsPrescribed;
            }

            else{
                //interpolate to find the value function at the foot of the characteristic
                nextVal = generalInterp(aU, aOC, aUYRG, aOCYRG, nextTimeIdx, dNext, vNext, tNext, aCurrentPos, aCurrentVel, aCurrentTime, currentControl);
            }

            //compute running cost integral
            double runningCostInt = runningCostIntegral(currentControl, aCurrentTime, gTau, gGamma, gC1, gC2, gC3);

            //for model 2, not in use
            /*if((gModel2) && (aProbOfLightChange == 1)){
                nextVal = delta;
            }*/

            currentValue = runningCostInt + (1 - aProbOfLightChange) * nextVal + aProbOfLightChange * delta;
            //currentValue = ((delta >= gInfty) && (gInUncertainGreenPhase)) ? gInfty : currentValue; //for model 2, not in use

        }

        //comparison with best value thus far
        if (currentValue < bestValueFnVal){
            assert(currentValue < gInfty);
            bestValueFnVal = currentValue;
            bestOCVal = currentControl;
            dNextOpt = dNext;
            vNextOpt = vNext;
            tNextOpt = tNext;
            assert(dNextOpt <= gDMax);
            assert(vNextOpt <= gVMax);
        }

    }//end loop over control

    //-----------------------------------------Golden Section Search----------------------------------------------
    if(gGoldenSection){
        double searchIntervalWidth = gControlSpacing;
        double aLeft = (bestOCVal - searchIntervalWidth < -gAlpha) ? -gAlpha : bestOCVal - searchIntervalWidth;
        if(gVelIndex == 0){
            aLeft = 0;
        }

        double aRight = (bestOCVal + searchIntervalWidth > gBeta) ? gBeta : bestOCVal + searchIntervalWidth;

        assert(aLeft >= -gAlpha);
        assert(aRight <= gBeta);

        struct sGSS gssValues = goldenSection(aU, aOC, aUYRG, aOCYRG, aUFirstPassYRG, aOCFirstPassYRG, aCurrentTimeIndex, aLeft, aRight, aCurrentTime, aCurrentPos, aCurrentVel, aTau);

        bestValueFnVal = gssValues.valueFn;
        bestOCVal = gssValues.oc;
        dNextOpt = gssValues.dNext;
        vNextOpt = gssValues.vNext;
        tNextOpt = gssValues.tNext;
        assert(bestValueFnVal < gInfty);
    }//end GSS section

    S.valueFn = bestValueFnVal;
    S.oc = bestOCVal;
    S.dNextOfficial = dNextOpt;
    S.vNextOfficial = vNextOpt;
    S.tNextOfficial = tNextOpt;

    return S;
}

/*==============================================================================
* Golden Section Search (GSS)
*============================================================================*/
struct sGSS goldenSection (multiarray *aValueFunction, multiarray *aOC, multiarray *aUYRG, multiarray *aOCYRG, multiarray *aUYRGFirstPass, multiarray *aOCYRGFirstPass, short int aTimeIndex, double aLeftEndPt, double aRightEndPt, double aCurrentTime, double aDistCurrent, double aVelCurrent, double aTau){
    /*
     * Function: goldenSection()
     *
     * Arguments: value function and OC multiarrays, current timeslice index, left search interval endpoint, right search interval endpoint, current physical time, current position (aDistCurrent), current velocity, SL timestep (tau)
     *
     * Purpose: Executes golden section search for minimum function value in [aLeftEndPt, aRightEndPt]. Returns control value
     * a that yields this min.
     *
     * Returns: Struct containing value function, optimal control, dNext, vNext computed with the optimal control, and tNext computed with tau
     */

    sGSS gssValues;

    double computedMin;
    double valueFunctionValue;
    double dNextOpt;
    double vNextOpt;
    double tNextOpt;
    //--------------------------------Initialize search paramaters-----------------------------------------------------------
    double conditionalP = conditionalProb(aCurrentTime + gDT);

    double tol = 1e-5;
    double tauGS = 0.5 * (sqrt(5) - 1);
    double intervalLength = fabs(aRightEndPt - aLeftEndPt);
    int refinementIteration = 0;

    double leftIntervalEndPt = aLeftEndPt;
    double rightIntervalEndPt = aRightEndPt;

    double valueFn1;
    double valueFn2;

    double x1; //control value update 1
    double x2; //control value update 2

    short int nextTimeIdx = aTimeIndex + 1;

    //-------------------------------------First point-----------------------------------------------------------------------
    x1 = rightIntervalEndPt - tauGS * intervalLength;

    struct sNextCoordinates nextCoords1 = tauRefinement(aDistCurrent, aVelCurrent, aCurrentTime, x1, aTau);
    double dNext1 = nextCoords1.nextPos;
    double vNext1 = nextCoords1.nextVel;
    double tNext1 = nextCoords1.nextTime;
    double tau1 = nextCoords1.tauValue;
    gTau = tau1;

    if(gInStationarySolve){
        //Assign interp weights
        //Compute dNext, vNext
        double dRescaled = (dNext1 - gDTarget) / gDPos;
        double vRescaled = vNext1 / gDVel;

        double dC = gPosIndex;
        double dF = dC - 1;
        double vF = gVelIndex;
        double vC = vF + 1;

        int dCIdx = int(dC);
        int dFIdx = dCIdx - 1;
        int vFIdx = int(vF);
        int vCIdx = vFIdx + 1;

        //adjustment if j == gvMax -- downward bias
        vC = (gVelIndex == gVNum) ? gVelIndex : vC;
        vF = (gVelIndex == gVNum) ? gVelIndex - 1 : vF;
        vCIdx = (gVelIndex == gVNum) ? gVelIndex : vCIdx;
        vFIdx = (gVelIndex == gVNum) ? gVelIndex - 1 : vFIdx;
        double alphaBottomRight = ((dRescaled - dF) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
        double alphaBottomLeft = ((dC - dRescaled) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
        double alphaTopLeft = ((dC - dRescaled) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));
        double alphaTopRight = ((dRescaled - dF) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));

        //Left values have converged - take from 0 since these columns have been written with the converged solution
        double uBottomLeft = (*aValueFunction)[dFIdx][vFIdx][0];
        double uTopLeft = (*aValueFunction)[dFIdx][vCIdx][0];

        double uBottomRight = (*aValueFunction)[dCIdx][vFIdx][1]; //current guess from before //TODO: Double check indices
        double uTopRight = (*aValueFunction)[dCIdx][vCIdx][1]; //previous iteration update at j+1 saved in 1 spot

        double runningCost = runningCostIntegral(x1, 0, tau1, gGamma, gC1, gC2, gC3); //Note: current time arg is irrelevant

        if (gModifiedGS){
            double gamma = (1 / (1 - alphaBottomRight));

            if((dNext1 > gDMax) || (vNext1 > gVMax) || (dNext1 < gDTarget) || (vNext1 < 0)){
                valueFn1 = gInfty;
            }
            else{
                valueFn1 = gamma * (runningCost + alphaTopLeft * uTopLeft + alphaTopRight * uTopRight + alphaBottomLeft * uBottomLeft);
            }
        }

        else{
            valueFn1 = runningCost + alphaTopLeft * uTopLeft + alphaBottomLeft * uBottomLeft + alphaTopRight * uTopRight + alphaBottomRight * uBottomRight;
        }
    }//end stationary solve if

    else if ((gInTracer) && (aTau < gDT)){
        //double delta = terminalCost(aValueFunction, aOC, aUYRG, aOCYRG, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
        double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
        double nextVal = dvtInterpolation(aValueFunction, aOC, aUYRG, aOCYRG, vNext1, dNext1, aDistCurrent, aVelCurrent, aCurrentTime, tNext1, x1);


        double runningCostInt = runningCostIntegral(x1, aCurrentTime, tau1, gGamma, gC1, gC2, gC3);

        valueFn1 = runningCostInt + (1 - conditionalP) * nextVal + conditionalP * delta;
    }


    else{
        if ((gGreenRedGreenProb)){
            gHasBoundaryOn = (tNext1 >= gTG) ? false : true;
        }

        bool allowedLocPt1 = isAllowed(vNext1, dNext1, tNext1);

        double runningCostInt = runningCostIntegral(x1, aCurrentTime, tau1, gGamma, gC1, gC2, gC3);
        double vFIsPrescribed1 = prescribedValueFunction(aValueFunction, aOC, dNext1, vNext1, tNext1, aTau, x1);

        if (!allowedLocPt1){
            valueFn1 = gInfty;
        }

        else if (vFIsPrescribed1 >= 0){
            double delta = terminalCost(aValueFunction, aOC, aUYRGFirstPass, aOCYRGFirstPass, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
            //double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
            valueFn1 = runningCostInt + (1 - conditionalP) * vFIsPrescribed1 + delta * conditionalP;
            assert(valueFn1 >= 0);
        }

        else{
            double delta = terminalCost(aValueFunction, aOC, aUYRGFirstPass, aOCYRGFirstPass, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
            //double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
            double nextVal = generalInterp(aValueFunction, aOC, aUYRG, aOCYRG, nextTimeIdx, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1);
            /*if((gModel2) && (conditionalP == 1)){
                nextVal = delta;
            }*/
            valueFn1 = runningCostInt + (1 - conditionalP) * nextVal + conditionalP * delta;
            valueFn1 = (valueFn1 >= gInfty) ? gInfty : valueFn1;

            //Debugging:
            if(delta < 0){
                delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
                cout << delta << "\n";
            }
            assert(delta >= 0);
            assert(valueFn1 >= 0);
        }
    }

    //------------------------------------------------------Second point-----------------------------------------------------------------------
    x2 = leftIntervalEndPt + tauGS * intervalLength;

    struct sNextCoordinates nextCoords2 = tauRefinement(aDistCurrent, aVelCurrent, aCurrentTime, x2, aTau);
    double dNext2 = nextCoords2.nextPos;
    double vNext2 = nextCoords2.nextVel;
    double tNext2 = nextCoords2.nextTime;
    double tau2 = nextCoords2.tauValue;
    gTau = tau2;

    if(gInStationarySolve){
        //Assign interp weights
        //Compute dNext, vNext
        double dRescaled = (dNext2 - gDTarget) / gDPos;
        double vRescaled = vNext2 / gDVel;


        double dC = gPosIndex;
        double dF = dC - 1;
        double vF = gVelIndex;
        double vC = vF + 1;

        int dCIdx = int(dC);
        int dFIdx = dCIdx - 1;
        int vFIdx = int(vF);
        int vCIdx = vFIdx + 1;

        //adjustment if j == gvMax -- downward bias
        vC = (gVelIndex == gVNum) ? gVelIndex : vC;
        vF = (gVelIndex == gVNum) ? gVelIndex - 1 : vF;
        vCIdx = (gVelIndex == gVNum) ? gVelIndex : vCIdx;
        vFIdx = (gVelIndex == gVNum) ? gVelIndex - 1 : vFIdx;
        double alphaBottomRight = ((dRescaled - dF) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
        double alphaBottomLeft = ((dC - dRescaled) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
        double alphaTopLeft = ((dC - dRescaled) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));
        double alphaTopRight = ((dRescaled - dF) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));

        //Left values have converged - take from 0 since these columns have been written with the converged solution
        double uBottomLeft = (*aValueFunction)[dFIdx][vFIdx][0];
        double uTopLeft = (*aValueFunction)[dFIdx][vCIdx][0];

        double uBottomRight = (*aValueFunction)[dCIdx][vFIdx][1]; //current guess from before //TODO: Double check indices
        double uTopRight = (*aValueFunction)[dCIdx][vCIdx][1]; //previous iteration update at j+1 saved in 1 spot

        double runningCost = runningCostIntegral(x2, 0, tau2, gGamma, gC1, gC2, gC3); //Note: current time arg is irrelevant

        if (gModifiedGS){
            double gamma = (1 / (1 - alphaBottomRight));

            if((dNext2 > gDMax) || (vNext2 > gVMax) || (dNext2 < gDTarget) || (vNext2 < 0)){
                valueFn2 = gInfty;
            }
            else{
                valueFn2 = gamma * (runningCost + alphaTopLeft * uTopLeft + alphaTopRight * uTopRight + alphaBottomLeft * uBottomLeft);
            }
        }

        else{
            valueFn2 = runningCost + alphaTopLeft * uTopLeft + alphaBottomLeft * uBottomLeft + alphaTopRight * uTopRight + alphaBottomRight * uBottomRight;
        }
    }//end stationary solve if

    else if ((gInTracer) && (aTau < gDT)){
        //double delta = terminalCost(aValueFunction, dNext2, vNext2, tNext2, gMu, gDTarget);
        double nextVal = dvtInterpolation(aValueFunction, aOC, aUYRG, aOCYRG, vNext2, dNext2, aDistCurrent, aVelCurrent, aCurrentTime, tNext2, x2);
        //double delta = terminalCost(aValueFunction, aOC, aUYRG, aOCYRG, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
        double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
        /*if((gModel2) && (conditionalP == 1)){
            nextVal = delta;
        }*/
        double runningCostInt = runningCostIntegral(x2, aCurrentTime, tau2, gGamma, gC1, gC2, gC3);

        valueFn2 = runningCostInt + (1 - conditionalP) * nextVal + conditionalP * delta;
    }


    else{
        if ((gGreenRedGreenProb)){
            gHasBoundaryOn = (tNext2 >= gTG) ? false : true;
        }

        bool allowedLocPt2 = isAllowed(vNext2, dNext2, tNext2);

        double runningCostInt2 = runningCostIntegral(x2, aCurrentTime, tau2, gGamma, gC1, gC2, gC3);
        double vFIsPrescribed2 = prescribedValueFunction(aValueFunction, aOC, dNext2, vNext2, tNext2, aTau, x2); //NOTE tau change!

        if (!allowedLocPt2){
            valueFn2 = gInfty;
        }

        else if (vFIsPrescribed2 >= 0){
            //double delta = terminalCost(aValueFunction, aOC, aUYRG, aOCYRG, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
            double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
            valueFn2 = runningCostInt2 + (1 - conditionalP) * vFIsPrescribed2 + delta * conditionalP;
            assert(valueFn2 >= 0);
        }

        else{
            //assert(gTau == gDT);
            //double delta = terminalCost(aValueFunction, dNext2, vNext2, tNext2, gMu, gDTarget);
            double delta = terminalCost(aValueFunction, aOC, aUYRGFirstPass, aOCYRGFirstPass, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
            //double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
            double nextVal = generalInterp(aValueFunction, aOC, aUYRG, aOCYRG, nextTimeIdx, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2);
            /*if((gModel2) && (conditionalP == 1)){
                nextVal = delta;
            }*/
            valueFn2 = runningCostInt2 + (1 - conditionalP) * nextVal + conditionalP * delta;
            valueFn2 = (valueFn2 >= gInfty) ? gInfty : valueFn2;
            assert(valueFn2 >= 0);
        }
    }

    //----------------------------While loop to find min-----------------------------------------------------------------------------------------

    while (intervalLength > tol){

        if(valueFn2 > valueFn1){
            rightIntervalEndPt = x2;
            x2 = x1;
            intervalLength = fabs(rightIntervalEndPt - leftIntervalEndPt);
            x1 = rightIntervalEndPt - tauGS * intervalLength;
            //Reset vf2
            valueFn2 = valueFn1;

            //compute new value function 1
            sNextCoordinates nextCoords1 = tauRefinement(aDistCurrent, aVelCurrent, aCurrentTime, x1, aTau);
            dNext1 = nextCoords1.nextPos;
            vNext1 = nextCoords1.nextVel;
            tNext1 = nextCoords1.nextTime;
            tau1 = nextCoords1.tauValue;
            gTau = tau1;

            if(gInStationarySolve){
                //Assign interp weights
                //Compute dNext, vNext
                double dRescaled = (dNext1 - gDTarget) / gDPos;
                double vRescaled = vNext1 / gDVel;


                double dC = gPosIndex;
                double dF = dC - 1;
                double vF = gVelIndex;
                double vC = vF + 1;

                int dCIdx = int(dC);
                int dFIdx = dCIdx - 1;
                int vFIdx = int(vF);
                int vCIdx = vFIdx + 1;

                //adjustment if j == gvMax -- downward bias
                vC = (gVelIndex == gVNum) ? gVelIndex : vC;
                vF = (gVelIndex == gVNum) ? gVelIndex - 1 : vF;
                vCIdx = (gVelIndex == gVNum) ? gVelIndex : vCIdx;
                vFIdx = (gVelIndex == gVNum) ? gVelIndex - 1 : vFIdx;
                double alphaBottomRight = ((dRescaled - dF) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
                double alphaBottomLeft = ((dC - dRescaled) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
                double alphaTopLeft = ((dC - dRescaled) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));
                double alphaTopRight = ((dRescaled - dF) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));

                //Left values have converged - take from 0 since these columns have been written with the converged solution
                double uBottomLeft = (*aValueFunction)[dFIdx][vFIdx][0];
                double uTopLeft = (*aValueFunction)[dFIdx][vCIdx][0];

                double uBottomRight = (*aValueFunction)[dCIdx][vFIdx][1]; //current guess from before //TODO: Double check indices
                double uTopRight = (*aValueFunction)[dCIdx][vCIdx][1]; //previous iteration update at j+1 saved in 1 spot

                double runningCost = runningCostIntegral(x1, 0, tau1, gGamma, gC1, gC2, gC3); //Note: current time arg is irrelevant

                if (gModifiedGS){
                    double gamma = (1 / (1 - alphaBottomRight));

                    if((dNext1 > gDMax) || (vNext1 > gVMax) || (dNext1 < gDTarget) || (vNext1 < 0)){
                        valueFn1 = gInfty;
                    }
                    else{
                        valueFn1 = gamma * (runningCost + alphaTopLeft * uTopLeft + alphaTopRight * uTopRight + alphaBottomLeft * uBottomLeft);
                    }
                }

                else{
                    valueFn1 = runningCost + alphaTopLeft * uTopLeft + alphaBottomLeft * uBottomLeft + alphaTopRight * uTopRight + alphaBottomRight * uBottomRight;
                }
            }//end stationary solve if


            else if ((gInTracer) && (aTau < gDT)){
                //double delta = terminalCost(aValueFunction, aOC, aUYRG, aOCYRG, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
                double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
                double nextVal = dvtInterpolation(aValueFunction, aOC, aUYRG, aOCYRG, vNext1, dNext1, aDistCurrent, aVelCurrent, aCurrentTime, tNext1, x1);

                /*if((gModel2) && (conditionalP == 1)){
                    nextVal = delta;
                }*/

                double runningCostInt = runningCostIntegral(x1, aCurrentTime, tau1, gGamma, gC1, gC2, gC3);

                valueFn1 = runningCostInt + (1 - conditionalP) * nextVal + conditionalP * delta;
            }


            else{
                if ((gGreenRedGreenProb)){
                    gHasBoundaryOn = (tNext1 >= gTG) ? false : true;
                }

                bool allowedLocPt1 = isAllowed(vNext1, dNext1, tNext1);

                double runningCostInt1 = runningCostIntegral(x1, aCurrentTime, tau1, gGamma, gC1, gC2, gC3);
                double vFIsPrescribed1 = prescribedValueFunction(aValueFunction, aOC, dNext1, vNext1, tNext1, aTau, x1);

                if (!allowedLocPt1){
                    valueFn1 = gInfty;
                }

                else if (vFIsPrescribed1 >= 0){
                    double delta = terminalCost(aValueFunction, aOC, aUYRGFirstPass, aOCYRGFirstPass, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
                    //double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
                    valueFn1 = runningCostInt1 + (1 - conditionalP) * vFIsPrescribed1 + delta * conditionalP;
                    assert(valueFn1 >= 0);
                }

                else{
                    assert(gTau == gDT);
                    double delta = terminalCost(aValueFunction, aOC, aUYRGFirstPass, aOCYRGFirstPass, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
                    //double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1, gMu, gDTarget);
                    double nextVal = generalInterp(aValueFunction, aOC, aUYRG, aOCYRG, nextTimeIdx, dNext1, vNext1, tNext1, aDistCurrent, aVelCurrent, aCurrentTime, x1);

                    /*if((gModel2) && (conditionalP == 1)){
                        nextVal = delta;
                    }*/

                    valueFn1 = runningCostInt1 + (1 - conditionalP) * nextVal + conditionalP * delta;
                    valueFn1 = (valueFn1 >= gInfty) ? gInfty : valueFn1;
                    assert(valueFn1 >= 0);
                }
            }

        }//end vf2 > vf1

        else{
            leftIntervalEndPt = x1;
            x1 = x2;
            intervalLength = fabs(rightIntervalEndPt - leftIntervalEndPt);
            x2 = leftIntervalEndPt + tauGS * intervalLength;
            valueFn1 = valueFn2;

            //compute new v(x2)
            sNextCoordinates nextCoords2 = tauRefinement(aDistCurrent, aVelCurrent, aCurrentTime, x2, aTau);
            dNext2 = nextCoords2.nextPos;
            vNext2 = nextCoords2.nextVel;
            tNext2 = nextCoords2.nextTime;
            tau2 = nextCoords2.tauValue;
            gTau = tau2;
            if(gInStationarySolve){
                //Assign interp weights
                //Compute dNext, vNext
                double dRescaled = (dNext2 - gDTarget) / gDPos;
                double vRescaled = vNext2 / gDVel;


                double dC = gPosIndex;
                double dF = dC - 1;
                double vF = gVelIndex;
                double vC = vF + 1;

                int dCIdx = int(dC);
                int dFIdx = dCIdx - 1;
                int vFIdx = int(vF);
                int vCIdx = vFIdx + 1;

                //adjustment if j == gvMax -- downward bias
                vC = (gVelIndex == gVNum) ? gVelIndex : vC;
                vF = (gVelIndex == gVNum) ? gVelIndex - 1 : vF;
                vCIdx = (gVelIndex == gVNum) ? gVelIndex : vCIdx;
                vFIdx = (gVelIndex == gVNum) ? gVelIndex - 1 : vFIdx;
                double alphaBottomRight = ((dRescaled - dF) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
                double alphaBottomLeft = ((dC - dRescaled) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
                double alphaTopLeft = ((dC - dRescaled) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));
                double alphaTopRight = ((dRescaled - dF) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));

                //Left values have converged - take from 0 since these columns have been written with the converged solution
                double uBottomLeft = (*aValueFunction)[dFIdx][vFIdx][0];
                double uTopLeft = (*aValueFunction)[dFIdx][vCIdx][0];

                double uBottomRight = (*aValueFunction)[dCIdx][vFIdx][1]; //current guess from before //TODO: Double check indices
                double uTopRight = (*aValueFunction)[dCIdx][vCIdx][1]; //previous iteration update at j+1 saved in 1 spot

                double runningCost = runningCostIntegral(x2, 0, tau2, gGamma, gC1, gC2, gC3); //Note: current time arg is irrelevant

                if (gModifiedGS){
                    double gamma = (1 / (1 - alphaBottomRight));

                    if((dNext2 > gDMax) || (vNext2 > gVMax) || (dNext2 < gDTarget) || (vNext2 < 0)){
                        valueFn2 = gInfty;
                    }
                    else{
                        valueFn2 = gamma * (runningCost + alphaTopLeft * uTopLeft + alphaTopRight * uTopRight + alphaBottomLeft * uBottomLeft);
                    }
                }

                else{
                    valueFn2 = runningCost + alphaTopLeft * uTopLeft + alphaBottomLeft * uBottomLeft + alphaTopRight * uTopRight + alphaBottomRight * uBottomRight;
                }
            }//end stationary solve if

            else if ((gInTracer) && (aTau < gDT)){
                double delta = terminalCost(aValueFunction, aOC, aUYRGFirstPass, aOCYRGFirstPass, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
                //double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
                double nextVal = dvtInterpolation(aValueFunction, aOC, aUYRG, aOCYRG, vNext2, dNext2, aDistCurrent, aVelCurrent, aCurrentTime, tNext2, x2);

                /*if((gModel2) && (conditionalP == 1)){
                    nextVal = delta;
                }*/

                double runningCostInt = runningCostIntegral(x2, aCurrentTime, tau2, gGamma, gC1, gC2, gC3);

                valueFn2 = runningCostInt + (1 - conditionalP) * nextVal + conditionalP * delta;
            }

            else{
                if ((gGreenRedGreenProb)){
                    gHasBoundaryOn = (tNext2 >= gTG) ? false : true;
                }

                bool allowedLocPt2 = isAllowed(vNext2, dNext2, tNext2);

                double runningCostInt2 = runningCostIntegral(x2, aCurrentTime, tau2, gGamma, gC1, gC2, gC3);
                double vFIsPrescribed2 = prescribedValueFunction(aValueFunction, aOC, dNext2, vNext2, tNext2, aTau, x2);

                if (!allowedLocPt2){
                    valueFn2 = gInfty;
                }

                else if (vFIsPrescribed2 >= 0){
                    double delta = terminalCost(aValueFunction, aOC, aUYRGFirstPass, aOCYRGFirstPass, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
                    //double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
                    valueFn2 = runningCostInt2 + (1 - conditionalP) * vFIsPrescribed2 + delta * conditionalP;
                    assert(valueFn2 >= 0);
                }

                else{
                    assert(gTau == gDT);
                    double delta = terminalCost(aValueFunction, aOC, aUYRGFirstPass, aOCYRGFirstPass, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
                    //double delta = terminalCost(aUYRG, aOCYRG, aUYRGFirstPass, aOCYRGFirstPass, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2, gMu, gDTarget);
                    double nextVal = generalInterp(aValueFunction, aOC, aUYRG, aOCYRG, nextTimeIdx, dNext2, vNext2, tNext2, aDistCurrent, aVelCurrent, aCurrentTime, x2);

                    /*if((gModel2) && (conditionalP == 1)){
                        nextVal = delta;
                    }*/

                    valueFn2 = runningCostInt2 + (1 - conditionalP) * nextVal + conditionalP * delta;
                    valueFn2 = (valueFn2 >= gInfty) ? gInfty : valueFn2;
                    assert(valueFn2 >= 0);
                }
            }
        }

        refinementIteration++;
        gNumParabPts = refinementIteration; //debugging, not relevant to any examples
    }

    //Set minimum on interval and continue
    if (valueFn2 > valueFn1){
        computedMin = x1;
        valueFunctionValue = valueFn1;
        dNextOpt = dNext1;
        vNextOpt = vNext1;
        tNextOpt = tNext1;

        vNextOpt = (fabs(vNextOpt) <= 1e-7) ? 0 : vNextOpt;
    }

    else {
        computedMin = x2;
        valueFunctionValue = valueFn2;
        dNextOpt = dNext2;
        vNextOpt = vNext2;
        tNextOpt = tNext2;

        vNextOpt = (fabs(vNextOpt) <= 1e-7) ? 0 : vNextOpt;

    }

    assert(valueFunctionValue >= 0);
    assert(valueFunctionValue < gInfty);
    assert(computedMin <= gBeta);
    gssValues.valueFn = valueFunctionValue;
    gssValues.oc = computedMin;
    gssValues.dNext = dNextOpt;
    gssValues.vNext = vNextOpt;
    gssValues.tNext = tNextOpt;
    return gssValues;
}

/*==============================================================================
* Determination if value function is prescribed
*============================================================================*/
double prescribedValueFunction(multiarray *aU, multiarray *aOC, double aD, double aV, double aT, double aTau, double aControl){
    /*
     Function: prescribedValueFunction()

     Arguments: value function and optimal control arrays, position (aD), velocity (aV), time (aT), SL timestep (aTau) and control value

     Purpose: Checks to see if the point is on a boundary or somewhere where value function is already known.
     If it is, return the proper value. If not, return negative value, indicating that interpolation needs to take place.

     Returns: Double representing the value function value at (aD, aV, aT)
     */

    double currentValue = -10;
    double targetDiff = fabs(aD - gDTarget);
    double tol = 1e-12;

    double timeToRedPhase = ((gTR - aT) > 0) ? (gTR - aT) : 0;
    double vCritical = gVMax - gBeta * timeToRedPhase; //velocity above which you'll reach vMax before TR if max accel

    //Quad form for velocity at intersection between max braking and max accel parabs
    double parabExtensionA = (gBeta + gAlpha) / (2 * gAlpha * gBeta); //a
    double parabExtensionB = -gVMax / gBeta; //b
    double parabExtensionC = pow(gVMax, 2) / (2 * gBeta) -timeToRedPhase * gVMax; //c
    double parabExtensionVelThresh = (-parabExtensionB + sqrt(pow(parabExtensionB, 2) - 4 *parabExtensionA * parabExtensionC)) / (2 * parabExtensionA); //quadratic formula

    double velocityThreshRoot = gAlpha * timeToRedPhase * (1 + sqrt(1 + (gBeta / gAlpha))); //velocity at intersection between max braking parab and max accel line
    double velocityThresh = (vCritical < velocityThreshRoot) ? parabExtensionVelThresh : velocityThreshRoot; //set velocity threshold
    velocityThresh = (timeToRedPhase > 0) ? velocityThresh : -tol;

    velocityThresh = (gFirstPass) ? -tol : velocityThresh;

    bool onParab = isOnParabPiece(aD, aV, aT);
    bool onMaxAccelBD = isOnMaxAccelBD(aD, aV, aT);

    double dParab = pow(aV, 2) / (2 * gAlpha);

    //tau refined on target
    if((aT >= gTG) || (!gGreenRedGreenProb)){
        if ((gTau < aTau) && (aD <= gDTarget)){
            assert(aD < 0);
            currentValue = targetBoundaryCost(aV, aT);
        }

        else if (targetDiff <= tol){
            currentValue = targetBoundaryCost(aV, aT);
        }
    }

    else{
        //tau refined, vehicle hits parabola
        if ((gTau < aTau) && (aT >= gTR) && (aT < gTG) && (aD >= dParab) && (gHasBoundaryOn)){
            currentValue = parabBoundaryCost(aU, aOC, aV, aT, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
        }//end tau refined

        //vehicle on max accel bd
        else if ((onMaxAccelBD) && (aV > velocityThresh) && (gHasBoundaryOn)){
            currentValue = maxAccelBoundaryCost(aU, aOC, aV, aT, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
        }

        //vehicle on parab bd
        else if ((onParab) && (aD != 0) && (aV != 0) && (gHasBoundaryOn)){
            currentValue = parabBoundaryCost(aU, aOC, aV, aT, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
        }

        //vehicle hits target
        else if ((gTau < aTau) && (aD <= gDTarget)){
            assert(aD < 0);
            currentValue = targetBoundaryCost(aV, aT);
        }

        //vehicle on target
        else if (targetDiff <= tol){
            currentValue = targetBoundaryCost(aV, aT);
        }

    }

    return currentValue;
}

/*==============================================================================
* Running Cost Integral
*============================================================================*/
double runningCostIntegral(double aControl, double aCurrentTime, double aTau, double aGamma, double aC1, double aC2, double aC3){
    /*
     Function: runningCostIntegral

     Arguments: control value, current physical time, SL timestep (tau), fuel cost parameter (aGamma, always set to 1), objective weights

     Purpose: Computes the integral of running cost. Also computes individual cost contributions when determining
     Pareto front.

     Returns: Double value for the integral of the running cost between t and t + aTau
     */

    double rCost = 0;
    double controlContribution = 0;
    controlContribution = (aControl > 0) ? aControl : 0; //[a]_+ term

    rCost = (aGamma * controlContribution * aC1 + aC2 * pow(aControl, 2) + aC3) * aTau; //integral of running cost between t and t+tau

    //Pareto = add up cost contributions
    if((gPareto) && (gSelectedTraj) && (gParetoCostComputation)){
        double fuelContribution = aGamma * controlContribution  * aTau;
        double discomfortContribution = pow(aControl, 2) * aTau;
        double timeContribution = aTau;

        gDiscomfortCost += pow(aControl, 2) * aTau;
        gFuelCost += aGamma * controlContribution  * aTau;
        gTimeCost += aTau;
        gTotalCost += (aC1 * fuelContribution) + (aC2 * discomfortContribution) + (aC3 * aTau);
    }

    return rCost;
}

/*==============================================================================
* Boundary cost along parabolic bd
*============================================================================*/
double parabBoundaryCost(multiarray *aU, multiarray *aOC, double aIntersectVel, double aIntersectT, double aTurningGreenT, double aGamma, double aC1, double aC2, double aC3, double aDeltaD, double aDeltaV, double aTargetPos){
    /*
     Function: parabBoundaryCost

     Arguments: value function and optimal control arrays, velocity on the boundary (aIntersectVel) time of intersection with boundary (aIntersectT),
     turning green time, fuel cost parameter (aGamma, always set to 1), objective weights, \Delta d, \Delta v, \Delta t, location of target

     Purpose: Computes the value function on the max-braking parabola d_alpha

     Returns: Double value function value
     */

    double bdCost = 0;

    //First, compute integral from bd intersect time to the time spent on parabola
    double stoppingTime = aIntersectVel / gAlpha; //time to come to a complete stop
    double timeRemaining = aTurningGreenT - aIntersectT; //time left until TG
    double tauParab = min(stoppingTime, timeRemaining); //time spent braking

    double runningCostInt = runningCostIntegral(-gAlpha, aIntersectT, tauParab, aGamma, aC1, aC2, aC3);

    //Next, determine the value function at dFinal, vFinal, Tg
    double dIntersect = pow(aIntersectVel, 2) / (2 * gAlpha);
    double dFinal = (tauParab == stoppingTime) ? 0 : dIntersect - (aIntersectVel * tauParab - 0.5 * gAlpha * pow(tauParab, 2));
    double vFinal = (tauParab == stoppingTime) ? 0 : aIntersectVel - gAlpha * tauParab;

    short int tIndex = short(int(round(aTurningGreenT / gDT)));
    //correction for when in terminal penalty
    if((gInTerminalPenalty) && (!gGreenRedGreenProb)){
        short int greenPhaseDurationIndices = short(int(round(gTYMax / gDT)));
        tIndex = tIndex - greenPhaseDurationIndices;
    }

    double terminalCost;

    if (tauParab == stoppingTime){
        double parkTime = aIntersectT + tauParab;
        double parkDuration = aTurningGreenT - parkTime;
        double parkedRC = runningCostIntegral(0, parkTime, parkDuration, aGamma, aC1, aC2, aC3);

        terminalCost = bilinearInterp(aU, aOC, tIndex, dFinal, 0, gTG, 0);

        bdCost = runningCostInt + parkedRC + terminalCost;
    }

    else{
        terminalCost = bilinearInterp(aU, aOC, tIndex, dFinal, vFinal, gTG, -gAlpha);
        bdCost = runningCostInt + terminalCost;
    }

    return bdCost;
}

/*==============================================================================
* Boundary cost along max accel line
*============================================================================*/
double maxAccelBoundaryCost(multiarray *aU, multiarray *aOC, double aIntersectVel, double aIntersectT, double aTurningRedT, double aGamma, double aC1, double aC2, double aC3, double aDeltaD, double aDeltaV, double aTargetPos){
    /*
     Function: maxAccelBoundaryCost

     Arguments: value function and optimal control arrays, velocity on the boundary (aIntersectVel) time of intersection with boundary (aIntersectT),
     turning green time, fuel cost parameter (aGamma, always set to 1), objective weights, \Delta d, \Delta v, \Delta t, location of target

     Purpose: Computes the value function on the max-acceleration portion of the boundary

     Returns: Double value function value on max accel boundary
     */

    double bdCost = 0;

    double terminalCost;

    //First, compute integral from bd intersect time to the time spent on max accel line
    //Determine dIntersect
    double dIntersect;

    double timeRemaining = aTurningRedT - aIntersectT;
    double maxVelTime = (gVMax - aIntersectVel) / gBeta;
    double timeAtMaxVel = aIntersectT + maxVelTime;
    double tauLinear = min(maxVelTime, timeRemaining);

    assert(aTurningRedT == gTR);
    assert(tauLinear >= 0);

    //debugging
    if (aIntersectVel == -0){
        aIntersectVel = 0;
    }

    if (timeRemaining > maxVelTime){
        dIntersect =  ((-pow(aIntersectVel, 2)) / (2 * gBeta)) + ((aIntersectVel * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeRemaining * gVMax);
    }

    else{
        dIntersect = (aIntersectVel * timeRemaining) + (0.5 * gBeta * pow(timeRemaining, 2));
    }

    double runningCostInt = runningCostIntegral(gBeta, aIntersectT, tauLinear, aGamma, aC1, aC2, aC3);

    //Next, determine the value function at dFinal, vFinal, Tg
    double vFinal = aIntersectVel + gBeta * tauLinear;

    assert(vFinal >= 0);

    short int tIndex = short(int(round(aTurningRedT / gDT))); //turning red index

    double turningRedTTerminalCost = aTurningRedT - gTYMax;
    if(gInTerminalPenalty){
        tIndex = short(int(round(turningRedTTerminalCost / gDT)));
    }

    double dFinal = dIntersect - ((aIntersectVel * tauLinear + 0.5 * gBeta * pow(tauLinear, 2)) + gVMax * (timeRemaining - tauLinear));

    if (dFinal != 0){
        double dFinalDiff = fabs(dFinal - 0);
        dFinal = (dFinalDiff <= 1e-13) ? 0 : dFinal;
    }
    assert(dFinal == 0);

    vFinal = aIntersectVel + gBeta * tauLinear;

    if (fabs(vFinal - gVMax) <= 1e-12){
        vFinal = gVMax;
    }

    terminalCost = bilinearInterp(aU, aOC, tIndex, dFinal, vFinal, gTR, gBeta);

    if (tauLinear == maxVelTime){
        terminalCost = bilinearInterp(aU, aOC, tIndex, dFinal, vFinal, gTR, 0);
        double coastDuration = timeRemaining - maxVelTime;
        double coastingRC = runningCostIntegral(0, coastDuration, coastDuration, aGamma, aC1, aC2, aC3);

        bdCost = runningCostInt + coastingRC + terminalCost;

    }

    else{
        bdCost = runningCostInt + terminalCost;
    }

    return bdCost;
}

/*==============================================================================
* Boundary cost along target at d = gDTarget
*============================================================================*/
double targetBoundaryCost(double aV, double aT){
    /*
     Function: targetBoundaryCost

     Arguments: velocity and time

     Purpose: Computes the cost at the target. Currently, we assume there's
     no additional cost incurred for reaching the target.

     Returns: Double value function value
     */

    double bdCost = 0;

    return bdCost;
}

/*==============================================================================
* Terminal cost
*============================================================================*/
double terminalCost (multiarray *aU, multiarray *aOC, multiarray *aUParab, multiarray *aOCParab, double aD, double aV, double aT, double aCurrentPos, double aCurrentVel, double aCurrentTime, double aControl, double aMu, double aTargetPos){
    /*
     Function: terminalCost

     Arguments: value function and optimal control arrays (including results from parab-only solve), next position (aD), next velocity (aV), next time (aT), current position, velocity, and time, control value, terminal penalty parameter (aMu - currently set to 1 and not relevant), location of target

     Purpose: Computes cost of termination at time aT. Only relevant during the uncertain green phase.
     */
    double termCost = 0;
    if (gUncertainGreenProb){

        //check if dNext, vNext out of bounds with control a
       if ((aD < gDTarget) || (aV > gVMax) || (aV < 0) || (aD > gDMax)){
            termCost = gInfty;
        }

       else{

           //Step 1: Figure out where cell is in relation to max accel line.
           double tymaxMinusDT = gTerminalT - gDT;

           /*gInTerminalPenalty = true;
           termCost = generalInterp(aU, aOC, aUParab, aOCParab, gNt, aD, aV, gTerminalT, aCurrentPos, aCurrentVel, tymaxMinusDT, aControl); //be very careful - 0 index is used to index the bottom slice on YRG but terminalT and terminalT-gDT is for the max accel line computation in uncert green
           gInTerminalPenalty = false;*/

           //Convert physical coordinates to grid coordinates
           double dRescaled = (aD - gDTarget) / gDPos;
           double vRescaled = aV / gDVel;

           //Grid coords
           short int dC = short(int(ceil(dRescaled)));
           short int dF = short(int(floor(dRescaled)));
           short int vC = short(int(ceil(vRescaled)));
           short int vF = short(int(floor(vRescaled)));

           //Bilinear interp
           //phi values at each corner of the interpolation cell
           double uA = (*aU)[dF][vC][gNt];
           double uB = (*aU)[dC][vC][gNt];
           double uC = (*aU)[dF][vF][gNt];
           double uD = (*aU)[dC][vF][gNt];

           //bilinear interp to get phi
           double beta = vRescaled - vF;
           double gamma = dRescaled - dF;

           //vertical interpolation
           double Q2 = (1-beta)*uC + beta*uA;
           double Q1 = (1-beta)*uD + beta*uB;

           termCost = gamma*Q1 + (1-gamma)*Q2;
       }
        gInTerminalPenalty = false;
    }

    else{
        double distFromTarget = aD - aTargetPos; //target somewhere to left of d = 0
        termCost = aMu * pow(distFromTarget, 2);
    }

    return termCost;
}

/*==============================================================================
 * Mollifying Terminal Condition - DEBUG, NOT IN USE - DELETE
 *============================================================================*/
void mollifyTerminalCondition (multiarray *aInitialValues, multiarray *aMollifiedSoln, short int aNumberOfTimesteps, double aDiffusionCoef){
    /*
     Solve u_t = gamma^2(u_dd + u_vv) with neumann BC's, subject to u(d,v,0) = aInitialValues(d,v,0)
     */
    cout << "In mollifier" << "\n";
    double deltaD = gDPos;
    double deltaV = gDVel;
    double deltaT = 1 / (2 * pow(aDiffusionCoef, 2) * ((1 / (pow(deltaD, 2))) + (1 / (pow(deltaV, 2))))); //TODO: Check just to make sure CFL condition satisfied. Will need to likely redo timestepping
    //double deltaT = 0.001;
    cout << "deltaT: " << deltaT << "\n";
    //Add one cell on each side, top and bottom
    short int totalNumPos = gDNum + 3;
    short int totalNumVel = gVNum + 3;

    //Multiarray declaration
    multiarray fullMollSoln(boost::extents[totalNumPos][totalNumVel][aNumberOfTimesteps + 1]); //array for values of value function at each point in domain

    //Note, 0 index is the left / lower ghost cell and totalNum - 1 index is the right / upper ghost cell
    //Prescribe initial values
    if(!gHorizDiffusionProbOnly){
        for(short int i = 0; i < totalNumPos; i++){     // distance loop
            for(short int j = 0; j < totalNumVel; j++){ // velocity loop
                if((i == 0) && (j >= 1) && (j <= gVNum + 1)){
                    fullMollSoln[i][j][0] = (*aInitialValues)[1][j - 1][0];
                }

                else if((i == totalNumPos - 1) && (j >= 1) && (j <= gVNum + 1)){
                    fullMollSoln[i][j][0] = (*aInitialValues)[gDNum - 1][j - 1][0];
                }

                else if ((j == 0) && (i >= 1) && (i <= gDNum + 1)){
                    fullMollSoln[i][j][0] = (*aInitialValues)[i - 1][1][0];
                }
                else if ((j == totalNumVel - 1) && (i >= 1) && (i <= gDNum + 1)){
                    fullMollSoln[i][j][0] = (*aInitialValues)[i - 1][gVNum - 1][0];
                }
            }
        }
        for(short int i = 1; i <= totalNumPos - 2; i++){     // distance loop
            for(short int j = 1; j <= totalNumVel - 2; j++){ // velocity loop
                fullMollSoln[i][j][0] = (*aInitialValues)[i-1][j-1][0];
            }
        }


        //Now, initialize top slices of red / yellow phase arrays
        for(short int k = 1; k <= aNumberOfTimesteps; k++){ // time loop
            for(short int i = 1; i <= totalNumPos - 2; i++){     // distance loop
                for(short int j = 1; j <= totalNumVel - 2; j++){ // velocity loop
                    double uijkp1;
                    double uip1jk = fullMollSoln[i+1][j][k-1];
                    double uijk = fullMollSoln[i][j][k-1];
                    double uim1jk = fullMollSoln[i-1][j][k-1];

                    double uijp1k = fullMollSoln[i][j+1][k-1];
                    double uijm1k = fullMollSoln[i][j-1][k-1];

                    double posDeriv = (uip1jk - 2 * uijk + uim1jk) / (pow(deltaD, 2));
                    double velDeriv = (uijp1k - 2 * uijk + uijm1k) / (pow(deltaV, 2));

                    uijkp1 = uijk + deltaT * pow(aDiffusionCoef, 2) * (posDeriv + velDeriv);

                    fullMollSoln[i][j][k] = uijkp1;

                }
            }

            //Now replace ghost points for periodic bc's
            for(short int i = 0; i <= totalNumPos - 1; i++){     // distance loop
                for(short int j = 0; j <= totalNumVel - 1; j++){ // velocity loop
                    if((i == 0) && (j >= 1) && (j <= gVNum + 1)){
                        fullMollSoln[i][j][k] = fullMollSoln[2][j][k];
                    }

                    else if((i == totalNumPos - 1) && (j >= 1) && (j <= gVNum + 1)){
                        fullMollSoln[i][j][k] = fullMollSoln[totalNumPos - 3][j][k];
                    }

                    else if ((j == 0) && (i >= 1) && (i <= gDNum + 1)){
                        fullMollSoln[i][j][k] = fullMollSoln[i][2][k];
                    }
                    else if ((j == totalNumVel - 1) && (i >= 1) && (i <= gDNum + 1)){
                        fullMollSoln[i][j][k] = fullMollSoln[i][totalNumVel - 3][k];
                    }
                }
            }
        }

        //Now, take interior of fullMollSoln and populate the aMollifiedSoln array with the final value
        for(short int i = 1; i <= totalNumPos - 2; i++){     // distance loop
            for(short int j = 1; j <= totalNumVel - 2; j++){ // velocity loop
                (*aMollifiedSoln)[i-1][j-1][0] = fullMollSoln[i][j][aNumberOfTimesteps];
            }
        }
    }

    else{
        multiarray fullMollSolnHoriz(boost::extents[totalNumPos][gVNum + 1][aNumberOfTimesteps + 1]);
        deltaT = pow(deltaD, 2) / (2 * pow(aDiffusionCoef, 2));

        //Set ghost pts
        for(short int i = 0; i < totalNumPos; i++){     // distance loop
            for(short int j = 0; j <= gVNum; j++){ // velocity loop
                if(i == 0){
                    fullMollSolnHoriz[i][j][0] = (*aInitialValues)[1][j][0];
                }

                else if(i == totalNumPos - 1){
                    fullMollSolnHoriz[i][j][0] = (*aInitialValues)[gDNum - 1][j][0];
                }

            }
        }

        for(short int i = 1; i <= totalNumPos - 2; i++){     // distance loop
            for(short int j = 0; j <= gVNum; j++){ // velocity loop
                fullMollSolnHoriz[i][j][0] = (*aInitialValues)[i-1][j][0];
            }
        }

        //Now, initialize top slices of red / yellow phase arrays
        for(short int k = 1; k <= aNumberOfTimesteps; k++){ // time loop
            for(short int j = 0; j <= gVNum; j++){     // distance loop
                for(short int i = 1; i <= totalNumPos - 2; i++){ // velocity loop
                    double uijkp1;
                    double uip1jk = fullMollSolnHoriz[i+1][j][k-1];
                    double uijk = fullMollSolnHoriz[i][j][k-1];
                    double uim1jk = fullMollSolnHoriz[i-1][j][k-1];

                    double posDeriv = (uip1jk - 2 * uijk + uim1jk) / (pow(deltaD, 2));

                    uijkp1 = uijk + deltaT * pow(aDiffusionCoef, 2) * (posDeriv);

                    fullMollSolnHoriz[i][j][k] = uijkp1;

                    if(i == 2){
                        fullMollSolnHoriz[0][j][k] = uijkp1;
                    }

                    else if (i == totalNumPos - 3){
                        fullMollSolnHoriz[totalNumPos - 1][j][k] = uijkp1;
                    }
                }
            }
        }

        for(short int i = 1; i <= totalNumPos - 2; i++){     // distance loop
            for(short int j = 0; j <= gVNum; j++){ // velocity loop
                (*aMollifiedSoln)[i-1][j][0] = fullMollSolnHoriz[i][j][aNumberOfTimesteps];
            }
        }

    }

}
/*==============================================================================
 * Infinity Norm
 *============================================================================*/
double infinityNorm(multiarray *aU, short int aPosIdx){
    /*
     Function: infinityNorm

     Arguments: value function and position index

     Purpose: Computes infinity norm at position slice aPosIdx for stationary GJ and GS solve.
     If aPosIdx < 0, computes infinity norm over entire domain.

     Returns: Double value of infinity norm.
     */

    double infNorm;

    short int iStartIdx = (aPosIdx < 0) ? 0 : aPosIdx;
    short int iEndIdx = (aPosIdx < 0) ? gDNum + 1 : iStartIdx + 1;

    double maxNorm = 0;
    for(short int i = iStartIdx; i < iEndIdx; i++) {     // distance loop
        for(short int j = 0; j < gVNum + 1; j++) { // velocity loop
            double err = fabs((*aU)[i][j][0] - (*aU)[i][j][1]);
            if(err > maxNorm){
                maxNorm = err;
            }
        }
    }

    infNorm = maxNorm;

    return  infNorm;
}

/*==============================================================================
 * Trajectory tracer
 *============================================================================*/
struct sOptimalTrajectory trajectoryTracer(multiarray *aValueFunction, multiarray *aOptimalControl, multiarray *aValueFunctionYRG, multiarray *aOptimalControlsYRG, multiarray *aValueFunctionYRGFirstPass, multiarray *aOptimalControlsYRGFirstPass, double aPosStart, double aVelStart, double aTimeStart, double aTau){
    /*
     Function: trajectoryTracer

     Arguments: value function and optimal control multiarrays, trajectory starting position, velocity, and time, and SL timestep (aTau)

     Purpose: Traces optimal trajectory starting from (aPosStart, aVelStart) at aTimeStart. Has option to trace the trajectory using
     value function results from the PDE solve or trace the trajectory by interpolating the optimal control values at each point in the domain
     determined during the PDE solve. Trajectory tracing ends once the vehicle has reached the target.

     NOTE: The default option is to trace the trajectory by interpolating the control values. This is computationally nicer and a bit cleaner
     overall. All examples in the manuscript rely on control interpolation to determine the optimal trajectories.

     Returns: A struct containing vectors storing the position, velocity, optimal control, and time values. When multiple possible
     turning yellow times are present, also returns vectors containing position, velocity, optimal control, and time for the
     trajectories traced associated with each possible turning yellow time.

     */

    /*-----------------------------------Vector initialization to hold stored values-----------------------------------------------------*/

    sOptimalTrajectory Traj;

    Traj.positionVals.push_back(aPosStart);
    Traj.velocityVals.push_back(aVelStart);
    Traj.timeVals.push_back(aTimeStart);

    if((gModel1Traj) && (gNumStops == 2)){
        Traj.s1Position.push_back(aPosStart);
        Traj.s1Velocity.push_back(aVelStart);
        Traj.s1Time.push_back(aTimeStart);
        Traj.s2Position.push_back(aPosStart);
        Traj.s2Velocity.push_back(aVelStart);
        Traj.s2Time.push_back(aTimeStart);
    }

    else if((gModel1Traj) && (gNumStops == 3)){
        Traj.s1Position.push_back(aPosStart);
        Traj.s1Velocity.push_back(aVelStart);
        Traj.s1Time.push_back(aTimeStart);

        Traj.s2Position.push_back(aPosStart);
        Traj.s2Velocity.push_back(aVelStart);
        Traj.s2Time.push_back(aTimeStart);

        Traj.s3Position.push_back(aPosStart);
        Traj.s3Velocity.push_back(aVelStart);
        Traj.s3Time.push_back(aTimeStart);

    }

    else if((gModel1Traj) && (gNumStops == 4)){
        Traj.s1Position.push_back(aPosStart);
        Traj.s1Velocity.push_back(aVelStart);
        Traj.s1Time.push_back(aTimeStart);

        Traj.s2Position.push_back(aPosStart);
        Traj.s2Velocity.push_back(aVelStart);
        Traj.s2Time.push_back(aTimeStart);

        Traj.s3Position.push_back(aPosStart);
        Traj.s3Velocity.push_back(aVelStart);
        Traj.s3Time.push_back(aTimeStart);

        Traj.s4Position.push_back(aPosStart);
        Traj.s4Velocity.push_back(aVelStart);
        Traj.s4Time.push_back(aTimeStart);
    }


    cout << "In trajectory tracer:" << "\n";

    gInTracer = true;

    int loopIndex = 0;

    bool hitTarget = false;
    bool loopDone = false;

    gHasBoundaryOn = false; //resetting red phase boolean from previous solves
    gTau = aTau;
    cout << "gTG:" << gTG << "\n";
    cout << "tau: " << gTau << "\n";

    double integratedCost = 0;
    double integratedRedPhaseCost = 0;
    double avgAccelRedPhase = 0;
    double redPhasePts = 0;

    /*--------------------------------------------------while loop to compute the trajectory------------------------------------------- */
    while(!loopDone){
        //current point
        double dCurrent = Traj.positionVals[loopIndex];
        //cout << "dCurrent: " << dCurrent << "\n";
        double vCurrent = Traj.velocityVals[loopIndex];
        //cout << "vCurrent: " << vCurrent <<"\n";
        double tCurrent = Traj.timeVals[loopIndex];
        //cout << "tCurrent: " << tCurrent << "\n";

        int timesliceIndex = round(tCurrent / gDT);
        double currentTCheck = timesliceIndex * gDT;
        double tCurrentDiff = fabs(currentTCheck - tCurrent);
        double dtTol = gDT / 1000;
        if(tCurrentDiff > dtTol){
            cout << "t current: " << tCurrent << " t current check: " << currentTCheck << "\n";
            cout << "timeslice idx: " << timesliceIndex << "\n";
        }

        assert(tCurrentDiff <= dtTol); //NOTE: off for tau < DT case!

        double dNextOfficial;
        double vNextOfficial;
        double tNextOfficial;

        double bestValueFnVal = gInfty;
        double bestOCVal = -gAccelInfty;

        double pStar = conditionalProb(tCurrent + gDT);
        if(pStar > 0){
            cout << "tcurrent for prob check: " << tCurrent << "\n";
        }

        //PW Boundary switches
        if(gUncertainGreenProb){
            if (tCurrent < gTYMax){
                gInUncertainGreenPhase = true;
                gHasBoundaryOn = false;
            }

            if ((tCurrent >= gTYMax)){
                gUncertainGreenProb = false;
                gGreenRedGreenProb = true;

                gHasBoundaryOn = (tCurrent < gTG) ? true : false;
            }
        }

        else if (gGreenRedGreenProb){
            gHasBoundaryOn = (tCurrent < gTG) ? true : false;
            assert(pStar == 0);
        }

        /* *********************************************Loop over controls ******************************************************** */
        bool startsInLegalRegion = isAllowed(vCurrent, dCurrent, tCurrent);
        double vfIsPrescribed = prescribedValueFunction(aValueFunction, aOptimalControl, dCurrent, vCurrent, tCurrent, aTau, 0);
        //Case 1: Point begins out of bounds, no need to go through control loop
        if (!startsInLegalRegion){
            bestValueFnVal = gInfty;
            bestOCVal = gAccelInfty;
            loopDone = true;
        }

        else if (vfIsPrescribed >= 0){
            cout << "hit prescribed vf" << "\n";
            cout << "loop idx: " << loopIndex << " timeslice idx: " <<timesliceIndex << "\n";
            bestValueFnVal = vfIsPrescribed;
            if(dCurrent <= gDTarget){
                bestOCVal = 0;
                dNextOfficial = gDTarget;
                hitTarget = true;
            }

            //glue on parabola when it is hit
            else{
                cout << "hit parabola" << "\n";
                cout << "tau: " << gTau << "\n";
                cout << "tCurrent: " << tCurrent << "\n";
                cout << "dCurrent: " << dCurrent << "\n";
                cout << "vCurrent: " << vCurrent << "\n";
                bool onParab = isOnParabPiece(dCurrent, vCurrent, tCurrent);
                //assert(onParab);
                double timeBrakingAlongParab = vCurrent / gAlpha;
                double timeSpentAtRest = gTG - timeBrakingAlongParab;
                double parabHitTime = tCurrent;
                bool endOfRedLightPhase = false;

                double dNextParab = Traj.positionVals[loopIndex];
                double vNextParab = Traj.velocityVals[loopIndex];
                double tNextParab = Traj.timeVals[loopIndex];

                //loop over and glue on
                while (!endOfRedLightPhase){
                    double dCurrentParab = Traj.positionVals[loopIndex];
                    double vCurrentParab = Traj.velocityVals[loopIndex];
                    double tCurrentParab = Traj.timeVals[loopIndex];

                    vNextParab = vCurrentParab - gAlpha * gDT;
                    dNextParab = pow(vNextParab, 2) / (2 * gAlpha);
                    tNextParab = tCurrentParab + gDT;
                    bestOCVal = -gAlpha;

                    if (tNextParab >= timeBrakingAlongParab + parabHitTime){
                        vNextParab = 0;
                        dNextParab = 0;
                        bestOCVal = 0;
                    }

                    dNextOfficial = dNextParab;
                    vNextOfficial = vNextParab;
                    tNextOfficial = tNextParab;

                    cout << "dNextOfficial: " << dNextOfficial << "vNextOff: " << vNextOfficial << "tNextOff: " << tNextOfficial << "bestOCVal: " << bestOCVal << "\n";

                    double integratedCostContribution = runningCostIntegral(bestOCVal, tCurrent, gDT, gGamma, gC1, gC2, gC3);

                    //save all but the last values here, last values saved and loop index advanced at bottom of loop
                    if (tNextOfficial < gTG){
                        Traj.positionVals.push_back(dNextOfficial);
                        Traj.velocityVals.push_back(vNextOfficial);
                        Traj.timeVals.push_back(tNextOfficial);
                        Traj.optimalControlVals.push_back(bestOCVal);
                        integratedCost += integratedCostContribution;

                        loopIndex++;
                    }

                    if(tNextParab >= gTG){
                        endOfRedLightPhase = true;
                    }
                }
            }
        }

        //Case 2: Point begins in bounds, can go through control loop
        else{
            if (gControlInterp){
                struct sOptimalValuesControlInterp optInterpVals;
                if(!gStationaryTraj){
                     optInterpVals = controlInterp(aValueFunction, aOptimalControl, aValueFunctionYRG, aOptimalControlsYRG, dCurrent, vCurrent, tCurrent, timesliceIndex, aTau, pStar);
                }

                else{
                    short int stationaryTimeIdx = gNt;
                    if(gInStationarySolve){
                        stationaryTimeIdx = 0;
                    }
                    optInterpVals = controlInterp(aValueFunction, aOptimalControl, aValueFunctionYRG, aOptimalControlsYRG, dCurrent, vCurrent, tCurrent, stationaryTimeIdx, aTau, pStar);

                }

                bestOCVal = optInterpVals.oc;
                dNextOfficial = optInterpVals.dNextOfficial;
                vNextOfficial = optInterpVals.vNextOfficial;
                tNextOfficial = optInterpVals.tNextOfficial;
                //cout << "tNextOfficial: " << tNextOfficial <<  " aTau: " << aTau << "\n";
            } //end control interp

            else{
                struct sOptimalValues optVals = controlLoop(aValueFunction, aOptimalControl, aValueFunctionYRG, aOptimalControlsYRG, aValueFunctionYRGFirstPass, aOptimalControlsYRGFirstPass, dCurrent, vCurrent, tCurrent, timesliceIndex, aTau, pStar);

                bestValueFnVal = optVals.valueFn;
                bestOCVal = optVals.oc;
                dNextOfficial = optVals.dNextOfficial;
                vNextOfficial = optVals.vNextOfficial;
                tNextOfficial = optVals.tNextOfficial;
                cout << "loop idx: " << loopIndex << "\n";
                assert(bestValueFnVal < gInfty);
            }
        }


        Traj.positionVals.push_back(dNextOfficial);
        Traj.velocityVals.push_back(vNextOfficial);
        Traj.timeVals.push_back(tNextOfficial);
        Traj.optimalControlVals.push_back(bestOCVal);

        if((gModel1Traj) && (gNumStops == 2)){
            if(tNextOfficial <= gTY1 + 1e-12){
                Traj.s1Position.push_back(dNextOfficial);
                Traj.s1Velocity.push_back(vNextOfficial);
                Traj.s1Time.push_back(tNextOfficial);
                Traj.s1Accel.push_back(bestOCVal);
            }

            Traj.s2Position.push_back(dNextOfficial);
            Traj.s2Velocity.push_back(vNextOfficial);
            Traj.s2Time.push_back(tNextOfficial);
            Traj.s2Accel.push_back(bestOCVal);
        }

        if((gModel1Traj) && (gNumStops == 3)){
            if(tNextOfficial <= gTY1 + 1e-12){
                Traj.s1Position.push_back(dNextOfficial);
                Traj.s1Velocity.push_back(vNextOfficial);
                Traj.s1Time.push_back(tNextOfficial);
                Traj.s1Accel.push_back(bestOCVal);
            }

            if(tNextOfficial <= gTY2 + 1e-12){
                Traj.s2Position.push_back(dNextOfficial);
                Traj.s2Velocity.push_back(vNextOfficial);
                Traj.s2Time.push_back(tNextOfficial);
                Traj.s2Accel.push_back(bestOCVal);
            }

            Traj.s3Position.push_back(dNextOfficial);
            Traj.s3Velocity.push_back(vNextOfficial);
            Traj.s3Time.push_back(tNextOfficial);
            Traj.s3Accel.push_back(bestOCVal);
        }

        if((gModel1Traj) && (gNumStops == 4)){
            if(tNextOfficial <= gTY1 + 1e-12){
                Traj.s1Position.push_back(dNextOfficial);
                Traj.s1Velocity.push_back(vNextOfficial);
                Traj.s1Time.push_back(tNextOfficial);
                Traj.s1Accel.push_back(bestOCVal);
            }

            if(tNextOfficial <= gTY2 + 1e-12){
                Traj.s2Position.push_back(dNextOfficial);
                Traj.s2Velocity.push_back(vNextOfficial);
                Traj.s2Time.push_back(tNextOfficial);
                Traj.s2Accel.push_back(bestOCVal);
            }

            if(tNextOfficial <= gTY3 + 1e-12){
                Traj.s3Position.push_back(dNextOfficial);
                Traj.s3Velocity.push_back(vNextOfficial);
                Traj.s3Time.push_back(tNextOfficial);
                Traj.s3Accel.push_back(bestOCVal);
            }

            Traj.s4Position.push_back(dNextOfficial);
            Traj.s4Velocity.push_back(vNextOfficial);
            Traj.s4Time.push_back(tNextOfficial);
            Traj.s4Accel.push_back(bestOCVal);
        }

        //pareto cost computation - compute cost of each objective:
        //Make sure you don't add contribution at t = gTG twice in Example 3 Pareto computation
        if(gPareto){
            gParetoCostComputation = true;
            double paretoCost = runningCostIntegral(bestOCVal, tCurrent, gTau, gGamma, gC1, gC2, gC3);
            gParetoCostComputation = false;
        }

        else{
            if((!gStationaryTraj)){
                gParetoCostComputation = true;
                double paretoCost = runningCostIntegral(bestOCVal, tCurrent, gTau, gGamma, gC1, gC2, gC3);
                gParetoCostComputation = false;
            }
            else if ((gStationaryTraj) && (tCurrent != 0)){
                gParetoCostComputation = true;
                double paretoCost = runningCostIntegral(bestOCVal, tCurrent, gTau, gGamma, gC1, gC2, gC3);
                gParetoCostComputation = false;
            }
        }

        double endingTime = (!gStationaryTraj) ? gTerminalT : 1000000;
        //Conditions to end loop
        if ((dNextOfficial <= gDTarget) || (hitTarget)|| (fabs(tNextOfficial - endingTime) <= 1e-12) || ((timesliceIndex + 1) > gNt)){

            loopDone = true;
            cout << "tNextOfficial" << tNextOfficial << "\n";
            cout << "dNextOfficial: " << dNextOfficial << "\n";
            cout << "oc: " << bestOCVal << "\n";
            if(((gModel1Traj) && (gNumStops == 2)) || (gModel2)){
                Traj.d0S2 = dNextOfficial;
                Traj.v0S2 = vNextOfficial;
            }

            if((gModel1Traj) && (gNumStops == 3)){
                Traj.d0S3 = dNextOfficial;
                Traj.v0S3 = vNextOfficial;
            }

            else if((gModel1Traj) && (gNumStops == 4)){
                Traj.d0S4 = dNextOfficial;
                Traj.v0S4 = vNextOfficial;
            }

        }

        loopIndex++; //advance loop index

        //Saving (d,v) at each of the Ti's in uncertain problem
        if((gModel1Traj) && (gNumStops == 2)){
            if (fabs(tNextOfficial - gTY1) <= 1e-12){
                Traj.d0S1 = dNextOfficial;
                Traj.v0S1 = vNextOfficial;
                cout << "tnext: " << tNextOfficial << "\n";
            }
        }

        else if((gModel1Traj) && (gNumStops == 3)){
            if (fabs(tNextOfficial - gTY1) <= 1e-12){
                Traj.d0S1 = dNextOfficial;
                Traj.v0S1 = vNextOfficial;
                cout << "tnext: " << tNextOfficial << "\n";
            }

            else if (fabs(tNextOfficial - gTY2) <= 1e-12){
                Traj.d0S2 = dNextOfficial;
                Traj.v0S2 = vNextOfficial;
                cout << "tnext: " << tNextOfficial << "\n";
            }
        }

        else if((gModel1Traj) && (gNumStops == 4)){
            if (fabs(tNextOfficial - gTY1) <= 1e-12){
                Traj.d0S1 = dNextOfficial;
                Traj.v0S1 = vNextOfficial;
                cout << "tnext: " << tNextOfficial << "\n";
            }

            else if (fabs(tNextOfficial - gTY2) <= 1e-12){
                Traj.d0S2 = dNextOfficial;
                Traj.v0S2 = vNextOfficial;
                cout << "tnext: " << tNextOfficial << "\n";
            }

            else if (fabs(tNextOfficial - gTY3) <= 1e-12){
                Traj.d0S3 = dNextOfficial;
                Traj.v0S3 = vNextOfficial;
                cout << "tnext: " << tNextOfficial << "\n";
            }
        }

        Traj.loopIdx = loopIndex;
        assert(loopIndex >=0);
    }//end while loop

    return Traj;

}

/*==============================================================================
 * All trajectories for Model 1, Two Lengths
 *============================================================================*/
void model1TwoLengthsTrajectories(multiarray *aValueFunctionUncertainGreen, multiarray *aOCUncertainGreen, multiarray *aYRGValueFunction, multiarray *aYRGOC, multiarray *aYRGFirstPassValueFunction, multiarray *aYRGFirstPassOC, multiarray *aStationaryValueFunction, multiarray *aStationaryOC, double aPosStart, double aVelStart, double aTimeStart, double aTau){

    /*
     Function: model1TwoLengthsTrajectories

     Arguments: value function and optimal control multiarrays, trajectory starting position, velocity, and time, and SL timestep (aTau)

     Purpose: Traces optimal trajectories for turning yellow time at T_1 and turning yellow time at T_2. Uses trajectoryTracer function
     and writes the optimal position, velocity, control and time for each to its own .txt file. See trajectoryTracer function for more
     information on the tracing itself.

     */

    struct sOptimalTrajectory optTraj = trajectoryTracer(aValueFunctionUncertainGreen, aOCUncertainGreen, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, aPosStart, aVelStart, aTimeStart, aTau);

    vector<double> optPos = optTraj.positionVals;
    vector<double> optVel = optTraj.velocityVals;
    vector<double> optTime = optTraj.timeVals;
    vector<double> optControls = optTraj.optimalControlVals;
    short int loopEndingIdx = optTraj.loopIdx;
    short int uncertainPhaseTrajEndIdx = loopEndingIdx;

    //Model 1 Trajectory Tracing - T_1
    vector<double> s1Pos = optTraj.s1Position;
    vector<double> s1Vel = optTraj.s1Velocity;
    vector<double> s1Time = optTraj.s1Time;
    vector<double> s1OptControls = optTraj.s1Accel;

    //Model 1 Trajectory Tracing - T_2
    vector<double> s2Pos = optTraj.s2Position;
    vector<double> s2Vel = optTraj.s2Velocity;
    vector<double> s2Time = optTraj.s2Time;
    vector<double> s2OptControls = optTraj.s2Accel;

    double s1InitialPos = optTraj.d0S1;
    double s1InitialVel = optTraj.v0S1;
    double s2InitialPos = optTraj.d0S2;
    double s2InitialVel = optTraj.v0S2;


    //Combine trajectory arrays
    if(gModel1Traj){
        vector<double> s1FullPosition;
        vector<double> s1FullVelocity;
        vector<double> s1FullOC;
        vector<double> s1FullTime;

        vector<double> s2FullPosition;
        vector<double> s2FullVelocity;
        vector<double> s2FullOC;
        vector<double> s2FullTime;

        gModel1Traj = false;
        gModel1TrajDetPortion = true;
        gUncertainGreenProb = false;
        gInUncertainGreenPhase = false;
        gGreenRedGreenProb = true;
        gModel1 = false;
        initGrid();
        gFirstTraj = true;

        cout << "Tracing S1" << "\n";
        cout << "d0: " << s1InitialPos << " v0: " << s1InitialVel << "\n";

        struct sOptimalTrajectory optTrajS1 = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, s1InitialPos, s1InitialVel, 0, aTau);

        vector<double> optPosDetS1 = optTrajS1.positionVals;
        vector<double> optVelDetS1 = optTrajS1.velocityVals;
        vector<double> optTimeDetS1 = optTrajS1.timeVals;
        vector<double> optControlsDetS1 = optTrajS1.optimalControlVals;

        gFirstTraj = false;
        gSelectedTraj = true;
        struct sOptimalTrajectory optTrajS2 = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, s2InitialPos, s2InitialVel, 0, aTau);
        gSelectedTraj = false;

        vector<double> optPosDetS2 = optTrajS2.positionVals;
        vector<double> optVelDetS2 = optTrajS2.velocityVals;
        vector<double> optTimeDetS2 = optTrajS2.timeVals;
        vector<double> optControlsDetS2 = optTrajS2.optimalControlVals;

        short int sizeS1Uncert = s1Time.size();
        short int sizeS1Det = optTimeDetS1.size();
        short int fullS1Size = sizeS1Det + sizeS1Uncert;

        short int sizeS2Uncert = s2Time.size();
        short int sizeS2Det = optTimeDetS2.size();
        short int fullS2Size = sizeS2Det + sizeS2Uncert;

        gModel1Traj = (gModel1) ? true : false;

        //combine both portions of the trajectory

        for(short int i = 0; i < sizeS1Uncert; i++) {
            double pos = s1Pos[i];
            double vel = s1Vel[i];
            double time = s1Time[i];
            double oc = s1OptControls[i];
            s1FullPosition.push_back(pos);
            s1FullVelocity.push_back(vel);
            s1FullTime.push_back(time);
            if(i < sizeS1Uncert - 1){
                s1FullOC.push_back(oc);
            }
        }

        //combine both portions of the trajectory
        for(short int i = 0; i < sizeS1Det; i++) {
            double pos = optPosDetS1[i];
            double vel = optVelDetS1[i];
            double time = optTimeDetS1[i] + gTY1;
            double oc = optControlsDetS1[i];
            assert(vel >= 0);
            //for stationary solve
            if(i < sizeS1Det - 1){
                s1FullOC.push_back(oc);
            }

            if(i > 0){
                s1FullPosition.push_back(pos);
                s1FullVelocity.push_back(vel);
                s1FullTime.push_back(time);
            }

        }
        //combine both portions of the trajectory
        for(short int i = 0; i < sizeS2Uncert; i++) {
            double pos = s2Pos[i];
            double vel = s2Vel[i];
            double time = s2Time[i];
            double oc = s2OptControls[i];
            s2FullPosition.push_back(pos);
            s2FullVelocity.push_back(vel);
            s2FullTime.push_back(time);
            if(i < sizeS2Uncert - 1){
                s2FullOC.push_back(oc);
            }

        }

        //combine both portions of the trajectory
        for(short int i = 0; i < sizeS2Det; i++){
            double pos = optPosDetS2[i];
            double vel = optVelDetS2[i];
            double time = optTimeDetS2[i] + gTY2;
            double oc = optControlsDetS2[i];

            if(i < sizeS2Det - 1){
                s2FullOC.push_back(oc);
            }

            if(i > 0){
                s2FullPosition.push_back(pos);
                s2FullVelocity.push_back(vel);
                s2FullTime.push_back(time);
            }
        }

        //If vehicle has not reached target by gTerminalT, trace the green light only portion of the traj to the target
        double traj1FinalTime = s1FullTime[fullS1Size - 2];
        double traj2FinalTime = s2FullTime[fullS2Size - 2];
        double traj1FinalPos = s1FullPosition[fullS1Size - 2];
        double traj2FinalPos = s2FullPosition[fullS2Size - 2];

        if((traj1FinalPos > gDTarget) && (traj1FinalTime >= gTerminalT)){
            double startingPos = s1FullPosition[fullS1Size - 2];
            double startingVel = s1FullVelocity[fullS1Size - 2];

            gStationaryTraj = true;
            gGreenRedGreenProb = false;
            gInStationarySolve = false;
            gUncertainGreenProb = false;

            struct sOptimalTrajectory optTrajStationary = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, startingPos, startingVel, 0, aTau);

            vector<double> optPosStat = optTrajStationary.positionVals;
            vector<double> optVelStat = optTrajStationary.velocityVals;
            vector<double> optTimeStat = optTrajStationary.timeVals;
            vector<double> optControlsStat = optTrajStationary.optimalControlVals;

            short int stationaryTrajSize = optTimeStat.size();

            for(short int i = 0; i < stationaryTrajSize; i++) {
                double pos = optPosStat[i];
                double vel = optVelStat[i];
                double time = optTimeStat[i] + gTerminalT + gTY1;
                double oc = optControlsStat[i];

                s1FullOC.push_back(oc);

                if(i > 0){
                    s1FullPosition.push_back(pos);
                    s1FullVelocity.push_back(vel);
                    s1FullTime.push_back(time);
                }
            }

            fullS1Size = fullS1Size + stationaryTrajSize;

        }

        if (traj2FinalPos > gDTarget){
            double startingPos = s2FullPosition[fullS2Size - 2];
            double startingVel = s2FullVelocity[fullS2Size - 2];
            gStationaryTraj = true;
            gGreenRedGreenProb = false;
            gInStationarySolve = false;
            gUncertainGreenProb = false;

            gSelectedTraj = true;
            struct sOptimalTrajectory optTrajStationary = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, startingPos, startingVel, 0, aTau);
            gSelectedTraj = false;

            vector<double> optPosStat = optTrajStationary.positionVals;
            vector<double> optVelStat = optTrajStationary.velocityVals;
            vector<double> optTimeStat = optTrajStationary.timeVals;
            vector<double> optControlsStat = optTrajStationary.optimalControlVals;

            short int stationaryTrajSize = optTimeStat.size();

            for(short int i = 0; i < stationaryTrajSize; i++) {
                double pos = optPosStat[i];
                double vel = optVelStat[i];
                double time = optTimeStat[i] + gTerminalT + gTY2;

                double oc = optControlsStat[i];

                s2FullOC.push_back(oc);

                if(i > 0){
                    s2FullPosition.push_back(pos);
                    s2FullVelocity.push_back(vel);
                    s2FullTime.push_back(time);
                }
            }

            fullS2Size = fullS2Size + stationaryTrajSize;
        }

        gModel1Traj = true; //back on for saving
        gFirstTraj = true;
        trajWriteToFile(s1FullPosition, s1FullVelocity, s1FullTime, s1FullOC, fullS1Size);
        cout << "written to file S1";
        gFirstTraj = false;
        trajWriteToFile(s2FullPosition, s2FullVelocity, s2FullTime, s2FullOC, fullS2Size);
        cout << "written to file s2";
    }
}

/*==============================================================================
 * All trajectories for Model 1, Three Lengths
 *============================================================================*/
void model1ThreeLengthsTrajectories(multiarray *aValueFunctionUncertainGreen, multiarray *aOCUncertainGreen, multiarray *aYRGValueFunction, multiarray *aYRGOC, multiarray *aYRGFirstPassValueFunction, multiarray *aYRGFirstPassOC, multiarray *aStationaryValueFunction, multiarray *aStationaryOC, double aPosStart, double aVelStart, double aTimeStart, double aTau){

    /*
     Function: model1ThreeLengthsTrajectories

     Arguments: value function and optimal control multiarrays, trajectory starting position, velocity, and time, and SL timestep (aTau)

     Purpose: Traces optimal trajectories for turning yellow time at T_1, turning yellow time at T_2, and turning yellow time at T_3.
     Uses trajectoryTracer function and writes the optimal position, velocity, control and time for each to its own .txt file.

     See trajectoryTracer function for more information on the tracing itself.

     */

    struct sOptimalTrajectory optTraj = trajectoryTracer(aValueFunctionUncertainGreen, aOCUncertainGreen, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, aPosStart, aVelStart, aTimeStart, aTau);

    vector<double> optPos = optTraj.positionVals;
    vector<double> optVel = optTraj.velocityVals;
    vector<double> optTime = optTraj.timeVals;
    vector<double> optControls = optTraj.optimalControlVals;
    short int loopEndingIdx = optTraj.loopIdx;
    short int uncertainPhaseTrajEndIdx = loopEndingIdx;


    //Model 1 Trajectory Tracing - T_1
    vector<double> s1Pos = optTraj.s1Position;
    vector<double> s1Vel = optTraj.s1Velocity;
    vector<double> s1Time = optTraj.s1Time;
    vector<double> s1OptControls = optTraj.s1Accel;

    //Model 1 Trajectory Tracing - T_2
    vector<double> s2Pos = optTraj.s2Position;
    vector<double> s2Vel = optTraj.s2Velocity;
    vector<double> s2Time = optTraj.s2Time;
    vector<double> s2OptControls = optTraj.s2Accel;

    //Model 1 Trajectory Tracing - T_3
    vector<double> s3Pos = optTraj.s3Position;
    vector<double> s3Vel = optTraj.s3Velocity;
    vector<double> s3Time = optTraj.s3Time;
    vector<double> s3OptControls = optTraj.s3Accel;

    double s1InitialPos = optTraj.d0S1;
    double s1InitialVel = optTraj.v0S1;

    double s2InitialPos = optTraj.d0S2;
    double s2InitialVel = optTraj.v0S2;

    double s3InitialPos = optTraj.d0S3;
    double s3InitialVel = optTraj.v0S3;

    cout << "Ending idx: " << loopEndingIdx << "\n";

    //Combine trajectory arrays
    if(gModel1Traj){
        vector<double> s1FullPosition;
        vector<double> s1FullVelocity;
        vector<double> s1FullOC;
        vector<double> s1FullTime;

        vector<double> s2FullPosition;
        vector<double> s2FullVelocity;
        vector<double> s2FullOC;
        vector<double> s2FullTime;

        vector<double> s3FullPosition;
        vector<double> s3FullVelocity;
        vector<double> s3FullOC;
        vector<double> s3FullTime;

        vector<double> s4FullPosition;
        vector<double> s4FullVelocity;
        vector<double> s4FullOC;
        vector<double> s4FullTime;

        gModel1Traj = false;
        gModel1TrajDetPortion = true;
        gUncertainGreenProb = false;
        gInUncertainGreenPhase = false;
        gGreenRedGreenProb = true;
        gModel1 = false;
        initGrid();
        gFirstTraj = true;

        cout << "Tracing S1" << "\n";
        cout << "d0: " << s1InitialPos << " v0: " << s1InitialVel << "\n";

        struct sOptimalTrajectory optTrajS1 = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, s1InitialPos, s1InitialVel, 0, aTau);

        vector<double> optPosDetS1 = optTrajS1.positionVals;
        vector<double> optVelDetS1 = optTrajS1.velocityVals;
        vector<double> optTimeDetS1 = optTrajS1.timeVals;
        vector<double> optControlsDetS1 = optTrajS1.optimalControlVals;

        gFirstTraj = false;

        //---------------Second traj-------------------------------------------------------
        gSecondTraj = true;
        gSelectedTraj = true;
        struct sOptimalTrajectory optTrajS2 = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, s2InitialPos, s2InitialVel, 0, aTau);
        gSelectedTraj = false;

        vector<double> optPosDetS2 = optTrajS2.positionVals;
        vector<double> optVelDetS2 = optTrajS2.velocityVals;
        vector<double> optTimeDetS2 = optTrajS2.timeVals;
        vector<double> optControlsDetS2 = optTrajS2.optimalControlVals;

        gSecondTraj = false;

        //---------------Third traj------------------------------------------------
        gThirdTraj = true;
        struct sOptimalTrajectory optTrajS3 = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, s3InitialPos, s3InitialVel, 0, aTau);
        gThirdTraj = false;
        vector<double> optPosDetS3 = optTrajS3.positionVals;
        vector<double> optVelDetS3 = optTrajS3.velocityVals;
        vector<double> optTimeDetS3 = optTrajS3.timeVals;
        vector<double> optControlsDetS3 = optTrajS3.optimalControlVals;

        //------------Piece together uncertain and det portions-----------------------------

        short int sizeS1Uncert = s1Time.size();
        short int sizeS1Det = optTimeDetS1.size();
        short int fullS1Size = sizeS1Det + sizeS1Uncert;

        short int sizeS2Uncert = s2Time.size();
        short int sizeS2Det = optTimeDetS2.size();
        short int fullS2Size = sizeS2Det + sizeS2Uncert;

        short int sizeS3Uncert = s3Time.size();
        short int sizeS3Det = optTimeDetS3.size();
        short int fullS3Size = sizeS3Det + sizeS3Uncert;

        gModel1Traj = (gModel1) ? true : false;

        //combine both portions of the trajectory
        cout << "size s1 uncert:" << sizeS1Uncert << "\n";
        for(short int i = 0; i < sizeS1Uncert; i++) {
            double pos = s1Pos[i];
            double vel = s1Vel[i];
            double time = s1Time[i];
            double oc = s1OptControls[i];
            s1FullPosition.push_back(pos);
            s1FullVelocity.push_back(vel);
            s1FullTime.push_back(time);
            if(i < sizeS1Uncert - 1){
                s1FullOC.push_back(oc);
            }
        }

        //combine both portions of the trajectory
        //for(short int i = 1; i < sizeS1Det; i++) {
        for(short int i = 0; i < sizeS1Det; i++) {
            double pos = optPosDetS1[i];
            double vel = optVelDetS1[i];
            double time = optTimeDetS1[i] + gTY1;
            double oc = optControlsDetS1[i];
            assert(vel >= 0);
            //for stationary solve

            if(i < sizeS1Det - 1){
                s1FullOC.push_back(oc);
            }

            if(i > 0){
                s1FullPosition.push_back(pos);
                s1FullVelocity.push_back(vel);
                s1FullTime.push_back(time);
            }

        }

        //Combine trajectory 2
        cout << "size s2 uncert:" << sizeS2Uncert << "\n";
        //combine both portions of the trajectory
        for(short int i = 0; i < sizeS2Uncert; i++) {
            double pos = s2Pos[i];
            double vel = s2Vel[i];
            double time = s2Time[i];
            double oc = s2OptControls[i];
            s2FullPosition.push_back(pos);
            s2FullVelocity.push_back(vel);
            s2FullTime.push_back(time);

            if(i < sizeS2Uncert - 1){
                s2FullOC.push_back(oc);
            }

        }

        //combine both portions of the trajectory
        for(short int i = 0; i < sizeS2Det; i++){
            double pos = optPosDetS2[i];
            double vel = optVelDetS2[i];
            double time = optTimeDetS2[i] + gTY2;
            double oc = optControlsDetS2[i];

            if(i < sizeS2Det - 1){
                s2FullOC.push_back(oc);
            }

            if(i > 0){
                s2FullPosition.push_back(pos);
                s2FullVelocity.push_back(vel);
                s2FullTime.push_back(time);
            }
        }

        //----------Combine trajectory 3----------------------------
        for(short int i = 0; i < sizeS3Uncert; i++) {
            double pos = s3Pos[i];
            double vel = s3Vel[i];
            double time = s3Time[i];
            double oc = s3OptControls[i];
            s3FullPosition.push_back(pos);
            s3FullVelocity.push_back(vel);
            s3FullTime.push_back(time);

            if(i < sizeS3Uncert - 1){
                s3FullOC.push_back(oc);
            }

        }

        //combine both portions of the trajectory
        for(short int i = 0; i < sizeS3Det; i++){
            double pos = optPosDetS3[i];
            double vel = optVelDetS3[i];
            double time = optTimeDetS3[i] + gTY3;
            double oc = optControlsDetS3[i];

            if(i < sizeS3Det - 1){
                s3FullOC.push_back(oc);
            }

            if(i > 0){
                s3FullPosition.push_back(pos);
                s3FullVelocity.push_back(vel);
                s3FullTime.push_back(time);
            }
        }

        //If vehicle has not reached target by gTerminalT, trace the green light only portion of the traj to the target
        double traj1FinalTime = s1FullTime[fullS1Size - 2];
        double traj2FinalTime = s2FullTime[fullS2Size - 2];
        double traj3FinalTime = s3FullTime[fullS3Size - 2];

        double traj1FinalPos = s1FullPosition[fullS1Size - 2];
        double traj2FinalPos = s2FullPosition[fullS2Size - 2];
        double traj3FinalPos = s3FullPosition[fullS3Size - 2];

        cout << "traj 2 final time: " << traj2FinalTime << "\n";

        if((traj1FinalPos > gDTarget) && (traj1FinalTime >= gTerminalT)){
            double startingPos = s1FullPosition[fullS1Size - 2];
            double startingVel = s1FullVelocity[fullS1Size - 2];
            cout << "stationary traj start pos: " << startingPos << " starting vel: " << startingVel << "\n";
            gStationaryTraj = true;
            initGrid();
            gGreenRedGreenProb = false;
            gInStationarySolve = false;
            gUncertainGreenProb = false;

            struct sOptimalTrajectory optTrajStationary = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, startingPos, startingVel, 0, aTau);

            vector<double> optPosStat = optTrajStationary.positionVals;
            vector<double> optVelStat = optTrajStationary.velocityVals;
            vector<double> optTimeStat = optTrajStationary.timeVals;
            vector<double> optControlsStat = optTrajStationary.optimalControlVals;

            short int stationaryTrajSize = optTimeStat.size();

            for(short int i = 0; i < stationaryTrajSize; i++) {
                double pos = optPosStat[i];
                double vel = optVelStat[i];
                double time = optTimeStat[i] + gTerminalT + gTY1;
                double oc = optControlsStat[i];

                s1FullOC.push_back(oc);

                if(i > 0){
                    s1FullPosition.push_back(pos);
                    s1FullVelocity.push_back(vel);
                    s1FullTime.push_back(time);
                }
            }

            fullS1Size = fullS1Size + stationaryTrajSize;

        }

        if(traj2FinalPos > gDTarget){
            double startingPos = s2FullPosition[fullS2Size - 2];
            double startingVel = s2FullVelocity[fullS2Size - 2];
            gStationaryTraj = true;
            initGrid();
            gGreenRedGreenProb = false;
            gInStationarySolve = false;
            gUncertainGreenProb = false;

            gSelectedTraj = true;

            struct sOptimalTrajectory optTrajStationary = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, startingPos, startingVel, 0, aTau);
            gSelectedTraj = false;

            vector<double> optPosStat = optTrajStationary.positionVals;
            vector<double> optVelStat = optTrajStationary.velocityVals;
            vector<double> optTimeStat = optTrajStationary.timeVals;
            vector<double> optControlsStat = optTrajStationary.optimalControlVals;

            short int stationaryTrajSize = optTimeStat.size();

            for(short int i = 0; i < stationaryTrajSize; i++) {// distance loop}
                double pos = optPosStat[i];
                double vel = optVelStat[i];
                double time = optTimeStat[i] + gTerminalT + gTY2;
                cout << "traj 2 time: " << time << "\n";
                double oc = optControlsStat[i];

                s2FullOC.push_back(oc);

                if(i > 0){
                    s2FullPosition.push_back(pos);
                    s2FullVelocity.push_back(vel);
                    s2FullTime.push_back(time);
                }

            }

            fullS2Size = fullS2Size + stationaryTrajSize;
        }

        if(traj3FinalPos > gDTarget){
            double startingPos = s3FullPosition[fullS3Size - 2];
            double startingVel = s3FullVelocity[fullS3Size - 2];
            gStationaryTraj = true;
            initGrid();
            gGreenRedGreenProb = false;
            gInStationarySolve = false;
            gUncertainGreenProb = false;

            gSelectedTraj = true;

            struct sOptimalTrajectory optTrajStationary = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, startingPos, startingVel, 0, aTau);
            gSelectedTraj = false;

            vector<double> optPosStat = optTrajStationary.positionVals;
            vector<double> optVelStat = optTrajStationary.velocityVals;
            vector<double> optTimeStat = optTrajStationary.timeVals;
            vector<double> optControlsStat = optTrajStationary.optimalControlVals;

            short int stationaryTrajSize = optTimeStat.size();

            for(short int i = 0; i < stationaryTrajSize; i++) {
                double pos = optPosStat[i];
                double vel = optVelStat[i];
                double time = optTimeStat[i] + gTerminalT + gTY3;
                cout << "traj 3 time: " << time << "\n";
                double oc = optControlsStat[i];

                s3FullOC.push_back(oc);

                if(i > 0){
                    s3FullPosition.push_back(pos);
                    s3FullVelocity.push_back(vel);
                    s3FullTime.push_back(time);
                }
            }

            fullS3Size = fullS3Size + stationaryTrajSize;
        }

        //----------------------Write output to a file-------------------------------------
        gModel1Traj = true;
        gFirstTraj = true;
        trajWriteToFile(s1FullPosition, s1FullVelocity, s1FullTime, s1FullOC, fullS1Size);
        cout << "written to file S1";
        gFirstTraj = false;

        gSecondTraj = true;
        trajWriteToFile(s2FullPosition, s2FullVelocity, s2FullTime, s2FullOC, fullS2Size);
        gSecondTraj = false;

        gThirdTraj = true;
        trajWriteToFile(s3FullPosition, s3FullVelocity, s3FullTime, s3FullOC, fullS3Size);
        gThirdTraj = false;

        cout << "written to file s2";
    }
}

/*==============================================================================
 * All trajectories for Model 1, Four Lengths
 *============================================================================*/
void model1FourLengthsTrajectories(multiarray *aValueFunctionUncertainGreen, multiarray *aOCUncertainGreen, multiarray *aYRGValueFunction, multiarray *aYRGOC, multiarray *aYRGFirstPassValueFunction, multiarray *aYRGFirstPassOC, multiarray *aStationaryValueFunction, multiarray *aStationaryOC, double aPosStart, double aVelStart, double aTimeStart, double aTau){

    /*
     Function: model1FourLengthsTrajectories

     Arguments: value function and optimal control multiarrays, trajectory starting position, velocity, and time, and SL timestep (aTau)

     Purpose: Traces optimal trajectories for turning yellow times T_1, T_2, T_3, and T_4. Uses trajectoryTracer function
     and writes the optimal position, velocity, control and time for each to its own .txt file. See trajectoryTracer function for more
     information on the tracing itself.

     */

    struct sOptimalTrajectory optTraj = trajectoryTracer(aValueFunctionUncertainGreen, aOCUncertainGreen, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, aPosStart, aVelStart, aTimeStart, aTau);

    vector<double> optPos = optTraj.positionVals;
    vector<double> optVel = optTraj.velocityVals;
    vector<double> optTime = optTraj.timeVals;
    vector<double> optControls = optTraj.optimalControlVals;
    short int loopEndingIdx = optTraj.loopIdx;
    short int uncertainPhaseTrajEndIdx = loopEndingIdx;

    trajWriteToFile(optPos, optVel, optTime, optControls, loopEndingIdx);

    //Model 1 Trajectory Tracing - T_1
    vector<double> s1Pos = optTraj.s1Position;
    vector<double> s1Vel = optTraj.s1Velocity;
    vector<double> s1Time = optTraj.s1Time;
    vector<double> s1OptControls = optTraj.s1Accel;

    //Model 1 Trajectory Tracing - T_2
    vector<double> s2Pos = optTraj.s2Position;
    vector<double> s2Vel = optTraj.s2Velocity;
    vector<double> s2Time = optTraj.s2Time;
    vector<double> s2OptControls = optTraj.s2Accel;

    //Model 1 Trajectory Tracing - T_3
    vector<double> s3Pos = optTraj.s3Position;
    vector<double> s3Vel = optTraj.s3Velocity;
    vector<double> s3Time = optTraj.s3Time;
    vector<double> s3OptControls = optTraj.s3Accel;

    //Model 1 Trajectory Tracing - T_4
    vector<double> s4Pos = optTraj.s4Position;
    vector<double> s4Vel = optTraj.s4Velocity;
    vector<double> s4Time = optTraj.s4Time;
    vector<double> s4OptControls = optTraj.s4Accel;

    double s1InitialPos = optTraj.d0S1;
    double s1InitialVel = optTraj.v0S1;

    double s2InitialPos = optTraj.d0S2;
    double s2InitialVel = optTraj.v0S2;

    double s3InitialPos = optTraj.d0S3;
    double s3InitialVel = optTraj.v0S3;

    double s4InitialPos = optTraj.d0S4;
    double s4InitialVel = optTraj.v0S4;


    cout << "Ending idx: " << loopEndingIdx << "\n";

    //Combine trajectory arrays
    if(gModel1Traj){
        vector<double> s1FullPosition;
        vector<double> s1FullVelocity;
        vector<double> s1FullOC;
        vector<double> s1FullTime;

        vector<double> s2FullPosition;
        vector<double> s2FullVelocity;
        vector<double> s2FullOC;
        vector<double> s2FullTime;

        vector<double> s3FullPosition;
        vector<double> s3FullVelocity;
        vector<double> s3FullOC;
        vector<double> s3FullTime;

        vector<double> s4FullPosition;
        vector<double> s4FullVelocity;
        vector<double> s4FullOC;
        vector<double> s4FullTime;

        gModel1Traj = false;
        gModel1TrajDetPortion = true;
        gUncertainGreenProb = false;
        gInUncertainGreenPhase = false;
        gGreenRedGreenProb = true;
        gModel1 = false;
        initGrid();
        gFirstTraj = true;

        cout << "Tracing S1" << "\n";
        cout << "d0: " << s1InitialPos << " v0: " << s1InitialVel << "\n";

        struct sOptimalTrajectory optTrajS1 = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, s1InitialPos, s1InitialVel, 0, aTau);

        vector<double> optPosDetS1 = optTrajS1.positionVals;
        vector<double> optVelDetS1 = optTrajS1.velocityVals;
        vector<double> optTimeDetS1 = optTrajS1.timeVals;
        vector<double> optControlsDetS1 = optTrajS1.optimalControlVals;

        gFirstTraj = false;

        //----------------Second traj-----------------------------------------------------
        gSecondTraj = true;
        gSelectedTraj = true;
        struct sOptimalTrajectory optTrajS2 = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, s2InitialPos, s2InitialVel, 0, aTau);
        gSelectedTraj = false;

        vector<double> optPosDetS2 = optTrajS2.positionVals;
        vector<double> optVelDetS2 = optTrajS2.velocityVals;
        vector<double> optTimeDetS2 = optTrajS2.timeVals;
        vector<double> optControlsDetS2 = optTrajS2.optimalControlVals;

        gSecondTraj = false;

        //------------Third traj-----------------------------------------------------
        gThirdTraj = true;
        struct sOptimalTrajectory optTrajS3 = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, s3InitialPos, s3InitialVel, 0, aTau);

        gThirdTraj = false;
        vector<double> optPosDetS3 = optTrajS3.positionVals;
        vector<double> optVelDetS3 = optTrajS3.velocityVals;
        vector<double> optTimeDetS3 = optTrajS3.timeVals;
        vector<double> optControlsDetS3 = optTrajS3.optimalControlVals;

        //-------------------Fourth traj---------------------------------------------
        gFourthTraj = true;
        struct sOptimalTrajectory optTrajS4 = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, s4InitialPos, s4InitialVel, 0, aTau);

        gFourthTraj = false;
        vector<double> optPosDetS4 = optTrajS4.positionVals;
        vector<double> optVelDetS4 = optTrajS4.velocityVals;
        vector<double> optTimeDetS4 = optTrajS4.timeVals;
        vector<double> optControlsDetS4 = optTrajS4.optimalControlVals;

        //---------------Piece together uncertain and det portions----------------------

        short int sizeS1Uncert = s1Time.size();
        short int sizeS1Det = optTimeDetS1.size();
        short int fullS1Size = sizeS1Det + sizeS1Uncert;

        short int sizeS2Uncert = s2Time.size();
        short int sizeS2Det = optTimeDetS2.size();
        short int fullS2Size = sizeS2Det + sizeS2Uncert;

        short int sizeS3Uncert = s3Time.size();
        short int sizeS3Det = optTimeDetS3.size();
        short int fullS3Size = sizeS3Det + sizeS3Uncert;

        short int sizeS4Uncert = s4Time.size();
        short int sizeS4Det = optTimeDetS4.size();
        short int fullS4Size = sizeS4Det + sizeS4Uncert;

        gModel1Traj = (gModel1) ? true : false;

        //combine both portions of the trajectory
        cout << "size s1 uncert:" << sizeS1Uncert << "\n";
        for(short int i = 0; i < sizeS1Uncert; i++) {
            double pos = s1Pos[i];
            double vel = s1Vel[i];
            double time = s1Time[i];
            double oc = s1OptControls[i];
            s1FullPosition.push_back(pos);
            s1FullVelocity.push_back(vel);
            s1FullTime.push_back(time);
            if(i < sizeS1Uncert - 1){
                s1FullOC.push_back(oc);
            }
        }

        //combine both portions of the trajectory
        for(short int i = 0; i < sizeS1Det; i++) {
            double pos = optPosDetS1[i];
            double vel = optVelDetS1[i];
            double time = optTimeDetS1[i] + gTY1;
            double oc = optControlsDetS1[i];
            assert(vel >= 0);

            if(i < sizeS1Det - 1){
                s1FullOC.push_back(oc);
            }

            if(i > 0){
                s1FullPosition.push_back(pos);
                s1FullVelocity.push_back(vel);
                s1FullTime.push_back(time);
            }

        }

        //Combine trajectory 2
        cout << "size s2 uncert:" << sizeS2Uncert << "\n";
        //combine both portions of the trajectory
        for(short int i = 0; i < sizeS2Uncert; i++) {
            double pos = s2Pos[i];
            double vel = s2Vel[i];
            double time = s2Time[i];
            double oc = s2OptControls[i];
            s2FullPosition.push_back(pos);
            s2FullVelocity.push_back(vel);
            s2FullTime.push_back(time);

            if(i < sizeS2Uncert - 1){
                s2FullOC.push_back(oc);
            }

        }

        //combine both portions of the trajectory
        for(short int i = 0; i < sizeS2Det; i++){
            double pos = optPosDetS2[i];
            double vel = optVelDetS2[i];
            double time = optTimeDetS2[i] + gTY2;
            double oc = optControlsDetS2[i];

            if(i < sizeS2Det - 1){
                s2FullOC.push_back(oc);
            }

            if(i > 0){
                s2FullPosition.push_back(pos);
                s2FullVelocity.push_back(vel);
                s2FullTime.push_back(time);
            }
        }

        //----------Combine trajectory 3----------------------------
        for(short int i = 0; i < sizeS3Uncert; i++) {
            double pos = s3Pos[i];
            double vel = s3Vel[i];
            double time = s3Time[i];
            double oc = s3OptControls[i];
            s3FullPosition.push_back(pos);
            s3FullVelocity.push_back(vel);
            s3FullTime.push_back(time);

            if(i < sizeS3Uncert - 1){
                s3FullOC.push_back(oc);
            }

        }

        //combine both portions of the trajectory
        for(short int i = 0; i < sizeS3Det; i++){
            double pos = optPosDetS3[i];
            double vel = optVelDetS3[i];
            double time = optTimeDetS3[i] + gTY3;
            double oc = optControlsDetS3[i];

            if(i < sizeS3Det - 1){
                s3FullOC.push_back(oc);
            }

            if(i > 0){
                s3FullPosition.push_back(pos);
                s3FullVelocity.push_back(vel);
                s3FullTime.push_back(time);
            }
        }

        //-----------------------------Combine fourth trajectory------------------------------------------------
        for(short int i = 0; i < sizeS4Uncert; i++) {
            double pos = s4Pos[i];
            double vel = s4Vel[i];
            double time = s4Time[i];
            double oc = s4OptControls[i];
            s4FullPosition.push_back(pos);
            s4FullVelocity.push_back(vel);
            s4FullTime.push_back(time);

            if(i < sizeS4Uncert - 1){
                s4FullOC.push_back(oc);
            }

        }

        //combine both portions of the trajectory
        for(short int i = 0; i < sizeS4Det; i++){
            double pos = optPosDetS4[i];
            double vel = optVelDetS4[i];
            double time = optTimeDetS4[i] + gTY4;
            double oc = optControlsDetS4[i];

            if(i < sizeS4Det - 1){
                s4FullOC.push_back(oc);
            }

            if(i > 0){
                s4FullPosition.push_back(pos);
                s4FullVelocity.push_back(vel);
                s4FullTime.push_back(time);
            }
        }

        //If vehicle has not reached target by gTerminalT, trace the green light only portion of the traj to the target
        double traj1FinalTime = s1FullTime[fullS1Size - 2];
        double traj2FinalTime = s2FullTime[fullS2Size - 2];
        double traj3FinalTime = s3FullTime[fullS3Size - 2];
        double traj4FinalTime = s4FullTime[fullS4Size - 2];

        double traj1FinalPos = s1FullPosition[fullS1Size - 2];
        double traj2FinalPos = s2FullPosition[fullS2Size - 2];
        double traj3FinalPos = s3FullPosition[fullS3Size - 2];
        double traj4FinalPos = s4FullPosition[fullS4Size - 2];

        cout << "traj 2 final time: " << traj2FinalTime << "\n";

        if((traj1FinalPos > gDTarget) && (traj1FinalTime >= gTerminalT)){
            double startingPos = s1FullPosition[fullS1Size - 2];
            double startingVel = s1FullVelocity[fullS1Size - 2];
            cout << "stationary traj start pos: " << startingPos << " starting vel: " << startingVel << "\n";
            gStationaryTraj = true;
            gGreenRedGreenProb = false;
            gInStationarySolve = false;
            gUncertainGreenProb = false;

            struct sOptimalTrajectory optTrajStationary = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, startingPos, startingVel, 0, aTau);

            vector<double> optPosStat = optTrajStationary.positionVals;
            vector<double> optVelStat = optTrajStationary.velocityVals;
            vector<double> optTimeStat = optTrajStationary.timeVals;
            vector<double> optControlsStat = optTrajStationary.optimalControlVals;

            short int stationaryTrajSize = optTimeStat.size();

            for(short int i = 0; i < stationaryTrajSize; i++) {
                double pos = optPosStat[i];
                double vel = optVelStat[i];
                double time = optTimeStat[i] + gTerminalT + gTY1;
                double oc = optControlsStat[i];

                s1FullOC.push_back(oc);

                if(i > 0){
                    s1FullPosition.push_back(pos);
                    s1FullVelocity.push_back(vel);
                    s1FullTime.push_back(time);
                }
            }

            fullS1Size = fullS1Size + stationaryTrajSize;

        }

        if(traj2FinalPos > gDTarget){
            double startingPos = s2FullPosition[fullS2Size - 2];
            double startingVel = s2FullVelocity[fullS2Size - 2];
            gStationaryTraj = true;
            gGreenRedGreenProb = false;
            gInStationarySolve = false;
            gUncertainGreenProb = false;

            gSelectedTraj = true;

            struct sOptimalTrajectory optTrajStationary = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, startingPos, startingVel, 0, aTau);
            gSelectedTraj = false;

            vector<double> optPosStat = optTrajStationary.positionVals;
            vector<double> optVelStat = optTrajStationary.velocityVals;
            vector<double> optTimeStat = optTrajStationary.timeVals;
            vector<double> optControlsStat = optTrajStationary.optimalControlVals;

            short int stationaryTrajSize = optTimeStat.size();

            for(short int i = 0; i < stationaryTrajSize; i++) {
                double pos = optPosStat[i];
                double vel = optVelStat[i];
                double time = optTimeStat[i] + gTerminalT + gTY2;
                cout << "traj 2 time: " << time << "\n";
                double oc = optControlsStat[i];

                s2FullOC.push_back(oc);

                if(i > 0){
                    s2FullPosition.push_back(pos);
                    s2FullVelocity.push_back(vel);
                    s2FullTime.push_back(time);
                }

            }

            fullS2Size = fullS2Size + stationaryTrajSize;
        }

        if(traj3FinalPos > gDTarget){
            double startingPos = s3FullPosition[fullS3Size - 2];
            double startingVel = s3FullVelocity[fullS3Size - 2];
            gStationaryTraj = true;
            gGreenRedGreenProb = false;
            gInStationarySolve = false;
            gUncertainGreenProb = false;

            gSelectedTraj = true;

            struct sOptimalTrajectory optTrajStationary = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, startingPos, startingVel, 0, aTau);
            gSelectedTraj = false;

            vector<double> optPosStat = optTrajStationary.positionVals;
            vector<double> optVelStat = optTrajStationary.velocityVals;
            vector<double> optTimeStat = optTrajStationary.timeVals;
            vector<double> optControlsStat = optTrajStationary.optimalControlVals;

            short int stationaryTrajSize = optTimeStat.size();

            for(short int i = 0; i < stationaryTrajSize; i++) {
                double pos = optPosStat[i];
                double vel = optVelStat[i];
                double time = optTimeStat[i] + gTerminalT + gTY3;
                cout << "traj 3 time: " << time << "\n";
                double oc = optControlsStat[i];

                s3FullOC.push_back(oc);

                if(i > 0){
                    s3FullPosition.push_back(pos);
                    s3FullVelocity.push_back(vel);
                    s3FullTime.push_back(time);
                }

            }

            fullS3Size = fullS3Size + stationaryTrajSize;
        }

        if(traj4FinalPos > gDTarget){
            double startingPos = s4FullPosition[fullS4Size - 2];
            double startingVel = s4FullVelocity[fullS4Size - 2];
            gStationaryTraj = true;
            gGreenRedGreenProb = false;
            gInStationarySolve = false;
            gUncertainGreenProb = false;

            gSelectedTraj = true;

            struct sOptimalTrajectory optTrajStationary = trajectoryTracer(aYRGValueFunction, aYRGOC, aYRGValueFunction, aYRGOC, aYRGFirstPassValueFunction, aYRGFirstPassOC, startingPos, startingVel, 0, aTau);
            gSelectedTraj = false;

            vector<double> optPosStat = optTrajStationary.positionVals;
            vector<double> optVelStat = optTrajStationary.velocityVals;
            vector<double> optTimeStat = optTrajStationary.timeVals;
            vector<double> optControlsStat = optTrajStationary.optimalControlVals;

            short int stationaryTrajSize = optTimeStat.size();

            for(short int i = 0; i < stationaryTrajSize; i++) {
                double pos = optPosStat[i];
                double vel = optVelStat[i];
                double time = optTimeStat[i] + gTerminalT + gTY4;
                cout << "traj 4 time: " << time << "\n";
                double oc = optControlsStat[i];

                s4FullOC.push_back(oc);

                if(i > 0){
                    s4FullPosition.push_back(pos);
                    s4FullVelocity.push_back(vel);
                    s4FullTime.push_back(time);
                }

            }

            fullS4Size = fullS4Size + stationaryTrajSize;
        }

        //------------------------------------------Write trajectories to file----------------------------------------------
        gModel1Traj = true; //back on for saving
        gFirstTraj = true;
        trajWriteToFile(s1FullPosition, s1FullVelocity, s1FullTime, s1FullOC, fullS1Size);
        cout << "written to file S1 \n";
        gFirstTraj = false;

        gSecondTraj = true;
        trajWriteToFile(s2FullPosition, s2FullVelocity, s2FullTime, s2FullOC, fullS2Size);
        cout << "written to file S2 \n";
        gSecondTraj = false;

        gThirdTraj = true;
        trajWriteToFile(s3FullPosition, s3FullVelocity, s3FullTime, s3FullOC, fullS3Size);
        cout << "written to file S3";
        gThirdTraj = false;

        gFourthTraj = true;
        trajWriteToFile(s4FullPosition, s4FullVelocity, s4FullTime, s4FullOC, fullS4Size);
        gFourthTraj = false;

        cout << "written to file S4 \n";
    }
}

/*==============================================================================
 * Pareto Grid
 *============================================================================*/
void paretoGrid(short int aNumPointsADim, short int aNumPointsBDim, double aC1, double aC2, double aC3, double aTrajStartPos, double aTrajStartVel, double aTrajStartTime){

    /*
     Function: paretoGrid()

     Arguments: number of points in A = C_3 dimension, number of points in B = C_1 / C_2 dimension, trajectory starting position, velocity, and time

     Purpose: Computes aNumPointsADim x aNumPointsBDim (A,B) points on the Pareto surface for (J_1 = total fuel cost, J_2 = total discomfort, J_3 = time to target). At each node, we store the computed (C_1, C_2, C_3) optimal costs for each of the three objectives. We are then able to use these values to determine the constrained Pareto front between two selected optimization objectives.

     Results are written to a .txt file for use in the j3ConstrainedParetoFront function.
     */

    gPareto = true;

    double selectedC3Val = 0.9;

    vector<double> discomfortCosts; //holds the discomfort costs, u - keLoss
    vector<double> fuelCosts; //holds final fuel costs
    vector<double> timeObj; //holds final times to target

    vector<double> c1WeightValues;
    vector<double> c2WeightValues;
    vector<double> c3WeightValues;

    multiarray paretoJ1Vals(boost::extents[aNumPointsADim+1][aNumPointsBDim+1][1]); //array for pareto J1 vals at each (A,B) node
    multiarray paretoJ2Vals(boost::extents[aNumPointsADim+1][aNumPointsBDim+1][1]); //array for pareto J2 vals at each (A,B) node
    multiarray paretoJ3Vals(boost::extents[aNumPointsADim+1][aNumPointsBDim+1][1]); //array for pareto J3 vals at each (A,B) node

    //Compute A values
    vector<double> nodeValuesA;

    double weightValASpacing = 1.0 / aNumPointsADim;

    for (short int i = 0; i < aNumPointsADim; i++){
        double nodeAVal = (i+1) * weightValASpacing;
        nodeValuesA.push_back(nodeAVal);

        cout << "A: " << nodeAVal << "\n";
    }

    //Compute B Values
    vector<double> nodeValuesB;
    double weightValBSpacing = 1.0 / (aNumPointsBDim - 1);
    for (short int i = 0; i < aNumPointsBDim; i++){
        double nodeBVal = i * weightValBSpacing;

        double squaredB = pow(nodeBVal, 2);
        nodeBVal = 1 - squaredB;

        nodeValuesB.push_back(nodeBVal);

        cout << "j: " << i << " B: " << nodeBVal << "\n";
    }

    //Pareto for each (A,B) pair
    for (short int i = 0; i < aNumPointsADim; i++){
        double nodeAVal = nodeValuesA[i];

        for (short int j = 0; j < aNumPointsBDim; j++){
            assert( j < aNumPointsBDim);
            double nodeBVal = nodeValuesB[j];

            double c3 = nodeAVal;
            double c1 = nodeBVal * (1 - c3);
            double c2 = 1 - (c1 + c3);

            c1WeightValues.push_back(c1);
            c2WeightValues.push_back(c2);
            c3WeightValues.push_back(c3);

            gC1 = c1;
            gC2 = c2;
            gC3 = c3;

            //cout << "Iteration: " << i << " c1: " << gC1 << " c2: " << gC2 << " c3: " << gC3 << "\n";

            //For uncertain, pick "hypothetical" turning yellow times
            //set aPosStart and aVelStart as global starting points for pareto
            if(gStationaryPareto){
                //cout << "stationary pareto" << "\n";
                greenLightOnlyProblem(aTrajStartPos, aTrajStartVel, aTrajStartTime);
            }

            else if (gYRGPareto){
                yellowRedProblem(aTrajStartPos, aTrajStartVel, aTrajStartTime);
            }

            else if(gUncertainPareto){
                yellowRedProblem(aTrajStartPos, aTrajStartVel, aTrajStartTime);
                uncertainGreenProblem(aTrajStartPos, aTrajStartVel, aTrajStartTime);
            }

            //Append costs to arrays
            if(selectedC3Val == nodeAVal){
                discomfortCosts.push_back(gDiscomfortCost);
                fuelCosts.push_back(gFuelCost);
                timeObj.push_back(gTimeCost);
            }
            //Display costs
            /*cout << "Total discomfort: " << gDiscomfortCost << "\n";
            cout << "Total fuel cost: " << gFuelCost << "\n";
            cout << "Total time: " << gTimeCost << "\n";*/

            //Save optimal costs for each objectives at node (A,B)
            paretoJ1Vals[i][j][0] = gFuelCost;
            paretoJ2Vals[i][j][0] = gDiscomfortCost;
            paretoJ3Vals[i][j][0] = gTimeCost;

            //--------------Resetting the global costs for each objective-----------------------
            gDiscomfortCost = 0;
            gFuelCost = 0;
            gTimeCost = 0;

        } //end loop over B values
    }//end loop over A values

    //Write to file
    paretoWriteToFile(&paretoJ1Vals, &paretoJ2Vals, &paretoJ3Vals, discomfortCosts, fuelCosts, timeObj, c1WeightValues, c2WeightValues, c3WeightValues, aNumPointsADim, aNumPointsBDim);
}

/*==============================================================================
 * J3 Constrained Pareto Front Generation
 *============================================================================*/
void j3ConstrainedParetoFront(short int aNumPointsADim, short int aNumPointsBDim, double aJ3ConstantValue, double aTrajStartPos, double aTrajStartVel, double aTrajStartTime){
    /*
     Function: j3ConstrainedParetoFront

     Arguments: number of points in A = C_3 dimension, number of points in B = C_1 / C_2 dimension, maximum time to target allowance (aJ3ConstantValue) trajectory starting position, velocity, and time

     Requires: j1_grid_stationary_row_maj.txt, j2_grid_stationary_row_maj.txt, j3_grid_stationary_row_maj.txt to read into multiarrays for use.

     Purpose: Takes in grid of (C_1, C_2, C_3) optimal costs for each of the three objectives at each (A,B) node. Uses the bisection method to determine the Pareto optimal J_1 and J_2 value for the indicated maximum time to target allowance.

     Results are written to a .txt file for plotting in the ParetoPlotting.ipynb notebook.
     */

    gDiscomfortCost = 0;
    gFuelCost = 0;
    gTimeCost = 0;

    gInConstrainedPareto = true;
    gPareto = true;

    //Init arrays
    multiarray paretoJ1Vals(boost::extents[aNumPointsADim+1][aNumPointsBDim+1][1]); //array for pareto J1 vals at each (A,B) node
    multiarray paretoJ2Vals(boost::extents[aNumPointsADim+1][aNumPointsBDim+1][1]); //array for pareto J2 vals at each (A,B) node
    multiarray paretoJ3Vals(boost::extents[aNumPointsADim+1][aNumPointsBDim+1][1]); //array for pareto J3 vals at each (A,B) node

    vector<double> discomfortCosts; //holds the discomfort costs
    vector<double> fuelCosts; //holds final fuel costs
    vector<double> timeObj; //holds final times to target

    vector<double> c1WeightValues;
    vector<double> c2WeightValues;
    vector<double> c3WeightValues;

    //---------------------------While loop over the refinement iterations---------------------------------------------------------
    //Compute A values
    vector<double> nodeValuesA;

    double weightValASpacing = 1.0 / aNumPointsADim;

    for (short int i = 0; i < aNumPointsADim; i++){
        double nodeAVal = (i+1) * weightValASpacing;
        nodeValuesA.push_back(nodeAVal);

        cout << "A: " << nodeAVal << "\n";
    }

    //Compute B Values
    vector<double> nodeValuesB;
    double weightValBSpacing = 1.0 / (aNumPointsBDim - 1);
    for (short int i = 0; i < aNumPointsBDim; i++){
        double nodeBVal = i * weightValBSpacing;

        double squaredB = pow(nodeBVal, 2);
        nodeBVal = 1 - squaredB;

        nodeValuesB.push_back(nodeBVal);

        cout << "j: " << i << " B: " << nodeBVal << "\n";
    }

    paretoReadFromFile(&paretoJ1Vals, &paretoJ2Vals, &paretoJ3Vals, aNumPointsADim, aNumPointsBDim);

    //Loop through B
    double tol = 1e-5;
    for (short int bValIdx = 0; bValIdx < aNumPointsBDim; bValIdx++){
        double bValue = nodeValuesB[bValIdx];

        double lowerAVal;
        double upperAVal;

        double j3AtLowerA;
        double j3AtUpperA;

        double j1AtLowerA;
        double j1AtUpperA;

        double j2AtLowerA;
        double j2AtUpperA;

        for(short int c3ValIdx = 1; c3ValIdx < aNumPointsADim; c3ValIdx++){
            double lowerAValCurrent = nodeValuesA[c3ValIdx - 1];
            double upperAValCurrent = nodeValuesA[c3ValIdx];

            double j3AtLowerACurrent = (paretoJ3Vals)[c3ValIdx - 1][bValIdx][0];
            double j3AtUpperACurrent = (paretoJ3Vals)[c3ValIdx][bValIdx][0];

            double j2AtLowerACurrent = (paretoJ2Vals)[c3ValIdx - 1][bValIdx][0];
            double j2AtUpperACurrent = (paretoJ2Vals)[c3ValIdx][bValIdx][0];

            double j1AtLowerACurrent = (paretoJ1Vals)[c3ValIdx - 1][bValIdx][0];
            double j1AtUpperACurrent = (paretoJ1Vals)[c3ValIdx][bValIdx][0];

            //As A increases, J3 decreases
            if ((aJ3ConstantValue <= j3AtLowerACurrent) && (aJ3ConstantValue >= j3AtUpperACurrent)){
                lowerAVal = lowerAValCurrent;
                upperAVal = upperAValCurrent;
                j3AtUpperA = j3AtUpperACurrent;
                j3AtLowerA = j3AtLowerACurrent;

                j2AtLowerA = j2AtLowerACurrent;
                j2AtUpperA = j2AtUpperACurrent;

                j1AtLowerA = j1AtLowerACurrent;
                j1AtUpperA = j1AtUpperACurrent;

                break;
            }
        }

        //-------------------------------------Bisection within interval-----------------------------------------
        double j3Diff = gInfty;
        double leftEndpoint = lowerAVal;
        double rightEndpoint = upperAVal;

        double leftJ3 = j3AtLowerA;
        double rightJ3 = j3AtUpperA;

        double savedJ1;
        double savedJ2;
        double savedJ3;
        double savedc3;
        double savedc2;
        double savedc1;

        bool bisectionLoopDone = false;

        short int iterationCounter = 0;

        while (!bisectionLoopDone){
            gDiscomfortCost = 0;
            gFuelCost = 0;
            gTimeCost = 0;

            double midpt = (leftEndpoint + rightEndpoint) * 0.5;

            //Compute HJB and costs with midpoint
            double currentC3 = midpt;
            double currentC1 = bValue * (1 - currentC3);
            double currentC2 = 1 - (currentC3 + currentC1);

            gC1 = currentC1;
            gC2 = currentC2;
            gC3 = currentC3;

            //cout << "Iteration: " << iterationCounter << " c1: " << gC1 << " c2: " << gC2 << " c3: " << gC3 << "\n";

            //Solve HJB and compute j1, j2, j3
            if(gStationaryPareto){
                greenLightOnlyProblem(aTrajStartPos, aTrajStartVel, aTrajStartTime);
            }

            else if (gYRGPareto){
                yellowRedProblem(aTrajStartPos, aTrajStartVel, aTrajStartTime);
            }

            else if(gUncertainPareto){
                yellowRedProblem(aTrajStartPos, aTrajStartVel, aTrajStartTime);
                uncertainGreenProblem(aTrajStartPos, aTrajStartVel, aTrajStartTime);
            }

            //Compare to the target fixed J3 value
            double j3DiffCurrent = fabs(gTimeCost - aJ3ConstantValue);
            double intervalHalfWidth = fabs(rightEndpoint - leftEndpoint)*0.5;
            if((j3DiffCurrent <= tol) || (intervalHalfWidth <= 1e-12)){
                bisectionLoopDone = true;
                savedc3 = midpt;
                savedJ1 = gFuelCost;
                savedJ2 = gDiscomfortCost;
                savedJ3 = gTimeCost;

                cout << "Total discomfort: " << gDiscomfortCost << "\n";
                cout << "Total fuel cost: " << gFuelCost << "\n";
                cout << "Total time: " << gTimeCost << "\n";

                savedc1 = currentC1;
                savedc2 = currentC2;
            }

            else{
                double midpointValDiff = gTimeCost - aJ3ConstantValue;
                double leftValDiff = leftJ3 - aJ3ConstantValue;
                double product = midpointValDiff * leftValDiff;
                if(product > 0){
                    //Set left val to c
                    leftEndpoint = midpt;
                    leftJ3 = gTimeCost;

                }
                else{
                    rightEndpoint = midpt;
                    rightJ3 = gTimeCost;
                }

            }//end value reset

            gDiscomfortCost = 0;
            gFuelCost = 0;
            gTimeCost = 0;

            iterationCounter++;

        }//end while

        //Save values
        discomfortCosts.push_back(savedJ2);
        fuelCosts.push_back(savedJ1);
        timeObj.push_back(savedJ3);

        //Save c3 value
        c3WeightValues.push_back(savedc3);
        c2WeightValues.push_back(savedc2);
        c1WeightValues.push_back(savedc1);
    }

    //write to file
    paretoWriteToFile(&paretoJ1Vals, &paretoJ2Vals, &paretoJ3Vals, discomfortCosts, fuelCosts, timeObj, c1WeightValues, c2WeightValues, c3WeightValues, aNumPointsADim, aNumPointsBDim);
}

/*==============================================================================
 * Write to file
 *============================================================================*/
void writeToFile (multiarray *aU, multiarray *aControls){
    /*
     Function: writeToFile()

     Purpose: Writes value function, optimal controls, and parameters used to a file.
     */

    std::ostringstream c2Obj;
    c2Obj << std::fixed;
    c2Obj << std::setprecision(2);
    c2Obj << gC2;
    std::string c2String = c2Obj.str();

    ofstream ufile;
    ofstream afile;
    ofstream vfile;

    if (gInStationarySolve){
        ufile.open("./output/stationary_valuefn.txt"); //value function values
        afile.open("./output/stationary_oc.txt"); //value function values
        vfile.open("./output/stationary_params.txt");
    }
    else if((!gUncertainGreenProb) && (!gInRGSolve)){
        ufile.open("./output/yrg_valuefn_c2_" + c2String + ".txt"); //value function values
        afile.open("./output/yrg_oc_c2_" + c2String + ".txt"); //optimal control values
        vfile.open("./output/yrg_params.txt");
    }

    else if((!gUncertainGreenProb) && (gInRGSolve)){
        ufile.open("./output/rg_valuefn.txt"); //value function values
        afile.open("./output/rg_oc.txt"); //value function values
        vfile.open("./output/rg_params.txt");
    }

    else{
        if(gNumStops == 2){
            double p1Val = gTYP1;
            double p2Val = gTYP2;

            std::ostringstream p1Obj;
            std::ostringstream p2Obj;

            p1Obj << std::fixed;
            p2Obj << std::fixed;

            p1Obj << std::setprecision(2);
            p2Obj << std::setprecision(2);

            p1Obj << p1Val;
            p2Obj << p2Val;

            std::string p1String = p1Obj.str();
            std::string p2String = p2Obj.str();

            ufile.open("./output/uncertain_green_valuefn_2L_c2_" + c2String + "_p1_" + p1String + "_p2_" + p2String + ".txt"); //value function values
            afile.open("./output/uncertain_green_oc_2L_c2_" + c2String + "_p1_" + p1String + "_p2_" + p2String + ".txt"); //optimal control values
            vfile.open("./output/uncertain_green_params_2L.txt");
        }
        else if (gNumStops == 3){
            double p1Val = gTYP1;
            double p2Val = gTYP2;
            double p3Val = gTYP3;

            std::ostringstream p1Obj;
            std::ostringstream p2Obj;
            std::ostringstream p3Obj;

            p1Obj << std::fixed;
            p2Obj << std::fixed;
            p3Obj << std::fixed;

            p1Obj << std::setprecision(2);
            p2Obj << std::setprecision(2);
            p3Obj << std::setprecision(2);

            p1Obj << p1Val;
            p2Obj << p2Val;
            p3Obj << p3Val;

            std::string p1String = p1Obj.str();
            std::string p2String = p2Obj.str();
            std::string p3String = p3Obj.str();

            ufile.open("./output/uncertain_green_valuefn_3L_c2_" + c2String + "_p1_" + p1String + "_p2_" + p2String + "_p3_" + p3String + ".txt"); //value function values
            afile.open("./output/uncertain_green_oc_3L_c2_" + c2String + "_p1_" + p1String + "_p2_" + p2String + "_p3_" + p3String + ".txt"); //optimal control values
            vfile.open("./output/uncertain_green_params_3L.txt");
        }
    }

    ofstream ubinfile;
    ubinfile.open("./output/green_light_binary.dat", ios::binary | ios::out);

    ofstream abinfile;
    abinfile.open("./output/green_light_binary.dat", ios::binary | ios::out);

  vfile << "gDPos" << " " << "gDVel" << " " << "gDT" << " " << "gDMax" << " " << "gVMax" << " " << "gDTarget" << " " << "gRedLightLength" << " " << "gDNum" << " " << "gVNum" << " " << "gNt" << " " << "gAlpha" << " " << "gBeta" << " " << "gTG" << " " << "gTR" << " " << "gTerminalT" << " " << "gAccelInfty" << " " << endl;
  vfile << gDPos << " " << gDVel << " " << gDT << " " << gDMax << " " << gVMax << " " << gDTarget << " " << gRedLightLength << " " << gDNum << " " << gVNum << " " << gNt << " " << gAlpha << " " << gBeta << " " << gTG << " " << gTR << " " << gTerminalT << " " << gAccelInfty << endl;
  vfile << "Probability Params:" << " " << "gTY1" << " " << "gTYMax" << " " << "gTYP1" << " " << "gTYP2" << " " << "gYellowDuration" << " " << " gRedDuration" << " " << "gTR (latest)" << " " << "gTG (latest)" << " "<< endl;
  vfile << gTY1 << " " << gTYMax << " " << gTYP1 << " " << gTYP2 << " " << gYellowDuration << " " << gRedDuration << " " << gTYMax + gYellowDuration << " " << gTYMax + gYellowDuration + gRedDuration << endl;
  vfile << "Sampling Params (for larger grids): " << "gSamplingTNumThresh" << " " << "gSliceSampleSpacing" << " "<< endl;
  vfile << gSamplingTNumThresh << " " << gSliceSampleSpacing << endl;


    vfile << gDPos << endl;
    vfile << gDVel << endl;
    vfile << gDT << endl;
    vfile << gDMax << endl;
    vfile << gVMax << endl;
    vfile << gDTarget << endl;
    vfile << gRedLightLength << endl;
    vfile << gDNum << endl;
    vfile << gVNum << endl;
    vfile << gNt << endl;
    vfile << gAlpha << endl;
    vfile << gBeta << endl;
    vfile << gTG << endl;
    vfile << gTR << endl;
    vfile << gSamplingTNumThresh << endl;
    vfile << gSliceSampleSpacing << endl;
    vfile.close();

   //ufile << gN << " " << gDMax << " " << gVMax << " " << gLightLength << " "<< gDNum << " " << gVNum <<endl;
   short int kEndIdx = (gInStationarySolve) ? 2 : gNt + 1;
   short int dEndIdx = (gInTSolve) ? (gDNum / 2) + 1 : gDNum + 1;
   for(short int k = gTimeEndIdx; k < kEndIdx; k++) {           // time
        for(short int i = 0; i < dEndIdx; i++) {      // distance
            for(short int j = 0; j < gVNum + 1; j++) {  // velocity

                //Sampling for realistic grid / larger number of timeslices. Write out every so many timeslices to file
                short int sliceMod = k % gSliceSampleSpacing;

                double val = (*aU)[i][j][k];

                if (gNt >= gSamplingTNumThresh){
                        if ((sliceMod == 0) || (k == 0)){

                        ufile << val;
                        ubinfile.write((char*)&val, sizeof(val));

                        double aval = (*aControls)[i][j][k];
                        afile << aval;
                        abinfile.write((char*)&aval, sizeof(aval));

                        if (j < gVNum) {
                            ufile << " ";
                            afile << " ";
                        }
                    }//end if write to file
                }//end sampling if

                else{
                    ufile << val;
                    ubinfile.write((char*)&val, sizeof(val));

                    double aval = (*aControls)[i][j][k];

                    afile << aval;
                    abinfile.write((char*)&aval, sizeof(aval));

                    if (j < gVNum) {
                        ufile << " ";
                        afile << " ";
                    }
                }//end else

            }
            ufile << endl;
            afile << endl;
        }
    }

    ufile.close();
    afile.close();
    ubinfile.close();
    abinfile.close();
}
/*=======================================================================================================================
 * Write to file - row major order
 *======================================================================================================================*/
void rowMajorWriteToFile(multiarray *aValueFunction, multiarray *aOptimalControls){
    /*
     Function: rowMajorWriteToFile

     Arguments: value function and optimal control multiarrays

     Purpose: Writes multiarray contents to file in row-major order.

     Order of terms: Stored in one long column by timeslice then by row.
     */

    std::ostringstream c2Obj;
    c2Obj << std::fixed;
    c2Obj << std::setprecision(2);
    c2Obj << gC2;
    std::string c2String = c2Obj.str();

    ofstream ufile;
    ofstream afile;

    if((gInStationarySolve) && (!gGridRef)){
        ufile.open("./output/stationary_value_fn_row_major_coarse.txt"); //value function values
        afile.open("./output/stationary_oc_row_major_coarse.txt"); //optimal control values
    }

    else if((gInStationarySolve) && (gGridRef) && (gGridRefFactor == 4)){
        ufile.open("./output/stationary_value_fn_row_major_fine.txt"); //value function values
        afile.open("./output/stationary_oc_row_major_fine.txt"); //optimal control values
    }

    else if((!gFirstPass) && (!gGridRef)){
        ufile.open("./output/value_fn_row_major_coarse_c2_" + c2String + ".txt"); //value function values
        afile.open("./output/oc_row_major_coarse_c2_" + c2String + ".txt"); //optimal control values
    }
    else if ((!gFirstPass) && (gGridRef) && (gGridRefFactor == 4)){
        ufile.open("./output/value_fn_row_major_fine_c2_" + c2String + ".txt"); //value function values
        afile.open("./output/oc_row_major_fine_c2_" + c2String + ".txt"); //optimal control values
    }
    else if ((gFirstPass) && (!gGridRef)){
        ufile.open("./output/value_fn_row_major_first_pass_coarse_c2_" + c2String + ".txt"); //value function values
        afile.open("./output/oc_row_major_first_pass_coarse_c2_" + c2String + ".txt"); //optimal control values
    }
    else if ((gFirstPass) && (gGridRef) && (gGridRefFactor == 4)){
        ufile.open("./output/value_fn_row_major_first_pass_fine_c2_" + c2String + ".txt"); //value function values
        afile.open("./output/oc_row_major_first_pass_fine_c2_" + c2String + ".txt"); //optimal control values
    }

    short int timeEndIdx = gNt;
    timeEndIdx = ((gInStationarySolve) || (gInConstrainedPareto)) ? 1 : gNt;

    for(short int k = 0; k < timeEndIdx + 1; k++) {           // time
         for(short int i = 0; i < gDNum + 1; i++) {      // distance
             for(short int j = 0; j < gVNum + 1; j++) {  // velocity

                 double val = (*aValueFunction)[i][j][k];
                 double aval = (*aOptimalControls)[i][j][k];

                 ufile << val << '\n';
                 afile << aval << '\n';

             }
         }
     }

     ufile.close();
     afile.close();
}

/*==============================================================================
 * HJB Solver read from file
 *============================================================================*/
void readFromFile(multiarray *aValueFunction, multiarray *aOptimalControls){
    /*
     Function: readFromFile

     Arguments: value function and optimal control multiarrays to hold data read in from row-major files.

     Note: Relies on the results being read-in from the row-major files

     Purpose: Reads the precomputed results for stage of interest from file to use as the terminal condition in uncertain problem
     */
    cout << "in read from file: " << "\n";
    std::ostringstream c2Obj;
    c2Obj << std::fixed;
    c2Obj << std::setprecision(2);
    c2Obj << gC2;
    std::string c2String = c2Obj.str();

    ifstream valueFnFile;
    ifstream ocFile;

    if((gInStationarySolve) && (!gGridRef)){
        valueFnFile.open("./output/stationary_value_fn_row_major_coarse.txt"); //value function values
        ocFile.open("./output/stationary_oc_row_major_coarse.txt"); //optimal control values
    }

    else if((gInStationarySolve) && (gGridRef) && (gGridRefFactor == 4)){
        valueFnFile.open("./output/stationary_value_fn_row_major_fine.txt"); //value function values
        ocFile.open("./output/stationary_oc_row_major_fine.txt"); //optimal control values
    }

    else if((!gFirstPass) && (!gGridRef)){
        cout << "corrected pass loading" << "\n";
        valueFnFile.open("./output/value_fn_row_major_coarse_c2_" + c2String + ".txt"); //value function values
        ocFile.open("./output/oc_row_major_coarse_c2_" + c2String + ".txt"); //optimal control values
    }
    else if ((!gFirstPass) && (gGridRef) && (gGridRefFactor == 4)){
        valueFnFile.open("./output/value_fn_row_major_fine_c2_" + c2String + ".txt"); //value function values
        ocFile.open("./output/oc_row_major_fine_c2_" + c2String + ".txt"); //optimal control values
    }
    else if ((gFirstPass) && (!gGridRef)){
        cout << "first pass loading" << "\n";
        valueFnFile.open("./output/value_fn_row_major_first_pass_coarse_c2_" + c2String + ".txt"); //value function values
        ocFile.open("./output/oc_row_major_first_pass_coarse_c2_" + c2String + ".txt"); //optimal control values
    }
    else if ((gFirstPass) && (gGridRef) && (gGridRefFactor == 4)){
        valueFnFile.open("./output/value_fn_row_major_first_pass_fine_c2_" + c2String + ".txt"); //value function values
        ocFile.open("./output/oc_row_major_first_pass_fine_c2_" + c2String + ".txt"); //optimal control values
    }

    short int timeEndIdx = gNt;
    timeEndIdx = (gInStationarySolve) ? 1 : gNt;

    for(short int k = 0; k < timeEndIdx + 1; k++) {           // time
         for(short int i = 0; i < gDNum + 1; i++) {      // distance
             for(short int j = 0; j < gVNum + 1; j++) {  // velocity
                 string ocVal;
                 string vfVal;
                 //cout << i << " " << j << " " << k << "\n";
                 getline(ocFile, ocVal, '\n');
                 getline(valueFnFile, vfVal, '\n');

                 //convert string to double
                 double ocDouble = stod(ocVal);
                 double vfDouble = stod(vfVal);

                 (*aValueFunction)[i][j][k] = vfDouble;
                 (*aOptimalControls)[i][j][k] = ocDouble;

             } //end vel loop
         } //end pos loop
    }//end time loop
    valueFnFile.close();
    ocFile.close();
}

/*==============================================================================
 * Write to file - trajectory
 *============================================================================*/
void trajWriteToFile (vector<double> aPosVec, vector<double> aVelVec, vector<double> aTimeVec, vector<double> aOCVals, short int aLoopEndIdx){
    /*
     Function: trajWriteToFile

     Arguments: vectors holding optimal position, velocity, optimal control values, and time, and the max-index
     value for loop indexing.

     Purpose: Writes results from trajectory tracer to .txt file in which
     the argument vectors are written out into lists.
     */

    //Check final time in trajectory to adjust loop end index accordingly
    short int finalIdx = aTimeVec.size() - 1;
    double finalTime = aTimeVec[finalIdx];

    /*-------------------------Write to file ------------------------------------------------------------------*/
    double d0Val = aPosVec[0];
    double v0Val = aVelVec[0];

    double p1Val = gTYP1;
    double p2Val = gTYP2;
    double p3Val = gTYP3;

    std::ostringstream d0StringObj;
    std::ostringstream v0StringObj;
    std::ostringstream p1Obj;
    std::ostringstream p2Obj;
    std::ostringstream p3Obj;
    std::ostringstream c2Obj;

    d0StringObj << std::fixed;
    v0StringObj << std::fixed;
    p1Obj << std::fixed;
    p2Obj << std::fixed;
    p3Obj << std::fixed;
    c2Obj << std::fixed;

    d0StringObj << std::setprecision(2);
    v0StringObj << std::setprecision(2);
    p1Obj << std::setprecision(2);
    p2Obj << std::setprecision(2);
    p3Obj << std::setprecision(2);
    c2Obj << std::setprecision(2);

    d0StringObj << d0Val;
    v0StringObj << v0Val;
    p1Obj << p1Val;
    p2Obj << p2Val;
    p3Obj << p3Val;
    c2Obj << gC2;

    std::string d0String = d0StringObj.str();
    std::string v0String = v0StringObj.str();
    std::string p1String = p1Obj.str();
    std::string p2String = p2Obj.str();
    std::string p3String = p3Obj.str();
    std::string c2String = c2Obj.str();

    ofstream tfile;

    if((gModel1Traj) && (gNumStops == 2)){
        if(gFirstTraj){
            cout << "OPEN S1, P1 PRIORITY!" << "P1: " << gTYP1 << "\n";
            tfile.open("./output/ug_c2_" + c2String + "_opt_trajectory_S1_d" + d0String + "v" + v0String + "_p1_" + p1String + "_p2_" + p2String + ".txt"); //optimal trajectory
        }
        else{
            tfile.open("./output/ug_c2_" + c2String + "_opt_trajectory_S2_d" + d0String + "v" + v0String + "_p1_" + p1String + "_p2_" + p2String + ".txt"); //optimal trajectory
        }
        aLoopEndIdx = aLoopEndIdx - 3;
    }

    else if((gModel1Traj) && (gNumStops == 3)){
        if(gFirstTraj){
            tfile.open("./output/ug_opt_trajectory_S1_d" + d0String + "v" + v0String + "_p1_" + p1String + "_p2_" + p2String + "_p3_" + p3String + ".txt"); //optimal trajectory
        }
        else if (gSecondTraj) {
            tfile.open("./output/ug_opt_trajectory_S2_d" + d0String + "v" + v0String + "_p1_" + p1String + "_p2_" + p2String + "_p3_" + p3String + ".txt"); //optimal trajectory
        }
        else if (gThirdTraj) {
            tfile.open("./output/ug_opt_trajectory_S3_d" + d0String + "v" + v0String + "_p1_" + p1String + "_p2_" + p2String + "_p3_" + p3String + ".txt"); //optimal trajectory
        }
        aLoopEndIdx = aLoopEndIdx - 3;
    }

    else if((gModel1Traj) && (gNumStops == 4)){
        if(gFirstTraj){
            tfile.open("./output/ug_opt_trajectory_S1_d" + d0String + "v" + v0String + ".txt"); //optimal trajectory
        }
        else if (gSecondTraj) {
            tfile.open("./output/ug_opt_trajectory_S2_d" + d0String + "v" + v0String + ".txt"); //optimal trajectory
        }
        else if (gThirdTraj) {
            tfile.open("./output/ug_opt_trajectory_S3_d" + d0String + "v" + v0String + ".txt"); //optimal trajectory
        }
        else if (gFourthTraj) {
            tfile.open("./output/ug_opt_trajectory_S4_d" + d0String + "v" + v0String + ".txt"); //optimal trajectory
        }
        aLoopEndIdx = finalIdx + 1;
    }

    else{
        if(gInStationarySolve){
            tfile.open("./output/stationary_opt_trajectory_d" + d0String + "v" + v0String + ".txt"); //optimal trajectory
        }

        else if ((gGreenRedGreenProb) && (gInRGSolve)){
            tfile.open("./output/rg_opt_trajectory_d" + d0String + "v" + v0String + ".txt"); //optimal trajectory

            aLoopEndIdx = finalIdx + 1;
        }

        else if ((gGreenRedGreenProb) && (!gInRGSolve)){
            tfile.open("./output/yrg_opt_trajectory_d" + d0String + "v" + v0String + ".txt"); //optimal trajectory
            aLoopEndIdx = finalIdx + 1;
        }

    }
    //--------------------------Optimal Trajectory write to file----------------------------------------------------

    for(int pos = 0; pos < aLoopEndIdx; pos++) {
        tfile <<  aPosVec[pos];
      if(pos != aLoopEndIdx - 1) {
        tfile << " ";
      }

    }
    tfile << endl;

    for(int vel = 0; vel < aLoopEndIdx; vel++) {
      tfile << aVelVec[vel];
        if(vel != aLoopEndIdx - 1) {
        tfile << " ";
      }
    }
    tfile << endl;

    for(int t = 0; t < aLoopEndIdx; t++) {
      tfile << aTimeVec[t];
      if(t != aLoopEndIdx - 1) {
        tfile << " ";
      }
    }
    tfile << endl;
    //tfile.close();
    for(int controlIdx = 0; controlIdx < aLoopEndIdx; controlIdx++) {
      tfile << aOCVals[controlIdx];
      if(controlIdx != aLoopEndIdx - 1) {
            tfile << " ";
      }
    }

    tfile.close();
}

/*==============================================================================
 * Write to file - pareto results
 *============================================================================*/
void paretoWriteToFile (multiarray *aJ1Array, multiarray *aJ2Array, multiarray *aJ3Array, vector<double> aDiscomfortCosts, vector<double> aFuelCosts, vector<double> aTimeCosts, vector<double> aC1Vals, vector<double> aC2Vals, vector<double> aC3Vals, short int aNumAPoints, short int aNumBPoints){
    /*
     Function: paretoWriteToFile

     Arguments: multiarrays containing the J1, J2, and J3 values at each node in the Pareto grid, vectors of optimal discomfort, fuel, and time to target from the J3 constrained Pareto front computation, number of points in A dimension, and number of points in B dimension.

     Purpose: Writes the vectors containing optimal discomfort, fuel, and time computed in the J3 constrained front function to file corresponding to the specified max time to target value (Example 1), and writes the J1, J2, and J3 grid values to a file in row major order.
     */

    ofstream tfile;
    ofstream j1File;
    ofstream j2File;
    ofstream j3File;
    ofstream j1FileRowMajor;
    ofstream j2FileRowMajor;
    ofstream j3FileRowMajor;

    if(gUncertainPareto){
        if (gJ3ConstrainedCurve1){
            tfile.open("./output/pareto_front_uncert_j3const_1.txt");
        }
        else if (gJ3ConstrainedCurve2){
            tfile.open("./output/pareto_front_uncert_j3const_2.txt");
        }
        else if (gJ3ConstrainedCurve3){
            tfile.open("./output/pareto_front_uncert_j3const_3.txt");
        }
        else{
            tfile.open("./output/pareto_front_uncert.txt");
        }
        j1File.open("./output/j1_grid_uncert.txt");
        j2File.open("./output/j2_grid_uncert.txt");
        j3File.open("./output/j3_grid_uncert.txt");

        if(!gInConstrainedPareto){
            j1FileRowMajor.open("./output/j1_grid_uncert_row_maj.txt");
            j2FileRowMajor.open("./output/j2_grid_uncert_row_maj.txt");
            j3FileRowMajor.open("./output/j3_grid_uncert_row_maj.txt");
        }
    }

    else if (gStationaryPareto){
        if (gJ3ConstrainedCurve1){
            tfile.open("./output/pareto_front_stationary_j3const_1.txt");
        }
        else if (gJ3ConstrainedCurve2){
            tfile.open("./output/pareto_front_stationary_j3const_2.txt");
        }
        else if (gJ3ConstrainedCurve3){
            tfile.open("./output/pareto_front_stationary_j3const_3.txt");
        }
        else{
            tfile.open("./output/pareto_front_stationary.txt");
        }

        j1File.open("./output/j1_grid_stationary.txt");
        j2File.open("./output/j2_grid_stationary.txt");
        j3File.open("./output/j3_grid_stationary.txt");

        if(!gInConstrainedPareto){
            j1FileRowMajor.open("./output/j1_grid_stationary_row_maj.txt");
            j2FileRowMajor.open("./output/j2_grid_stationary_row_maj.txt");
            j3FileRowMajor.open("./output/j3_grid_stationary_row_maj.txt");
        }
    }

    else if (gYRGPareto){
        if (gJ3ConstrainedCurve1){
            tfile.open("./output/pareto_front_yrg_j3const_1.txt");
        }
        else if (gJ3ConstrainedCurve2){
            tfile.open("./output/pareto_front_yrg_j3const_2.txt");
        }
        else if (gJ3ConstrainedCurve3){
            tfile.open("./output/pareto_front_yrg_j3const_3.txt");
        }
        else{
            tfile.open("./output/pareto_front_yrg.txt");
        }

        j1File.open("./output/j1_grid_yrg.txt");
        j2File.open("./output/j2_grid_yrg.txt");
        j3File.open("./output/j3_grid_yrg.txt");

        if(!gInConstrainedPareto){
            j1FileRowMajor.open("./output/j1_grid_yrg_row_maj.txt");
            j2FileRowMajor.open("./output/j2_grid_yrg_row_maj.txt");
            j3FileRowMajor.open("./output/j3_grid_yrg_row_maj.txt");
        }
    }

    short int numPts = aDiscomfortCosts.size();

    for(int i = 0; i < numPts; i++) {
        tfile <<  aDiscomfortCosts[i];
      if(i != numPts - 1) {
        tfile << " ";
      }
    }
    tfile << endl;

    for(int i = 0; i < numPts; i++) {
      tfile << aFuelCosts[i];
        if(i != numPts - 1) {
        tfile << " ";
      }
    }
    tfile << endl;

    for(int i = 0; i < numPts; i++) {
      tfile << aTimeCosts[i];
      if(i != numPts - 1) {
        tfile << " ";
      }
    }

    tfile << endl;

    for(int i = 0; i < numPts; i++) {
      tfile << aC1Vals[i];
        if(i != numPts - 1) {
        tfile << " ";
      }
    }
    tfile << endl;

    for(int i = 0; i < numPts; i++) {
      tfile << aC2Vals[i];
        if(i != numPts - 1) {
        tfile << " ";
      }
    }
    tfile << endl;

    for(int i = 0; i < numPts; i++) {
      tfile << aC3Vals[i];
        if(i != numPts - 1) {
        tfile << " ";
      }
    }

    tfile.close();

    //--------------------------Write Grid Values to File----------------------------------
    for(short int i = 0; i < aNumAPoints; i++) {      // distance
        for(short int j = 0; j < aNumBPoints; j++) {  // velocity

            double j1Val = (*aJ1Array)[i][j][0];
            double j2Val = (*aJ2Array)[i][j][0];
            double j3Val = (*aJ3Array)[i][j][0];

            j1File << j1Val;
            j2File << j2Val;
            j3File << j3Val;

            if (j < aNumBPoints - 1) {
                j1File << " ";
                j2File << " ";
                j3File << " ";

            }

        }
        j1File << endl;
        j2File << endl;
        j3File << endl;
    }

    j1File.close();
    j2File.close();
    j3File.close();
    /*------------------------------------------ROW MAJOR--------------------------------------------------------------*/

    if(!gInConstrainedPareto){
        for(short int i = 0; i < aNumAPoints; i++) { //A vals
            for(short int j = 0; j < aNumBPoints; j++) { //B vals

                double j1Val = (*aJ1Array)[i][j][0];
                double j2Val = (*aJ2Array)[i][j][0];
                double j3Val = (*aJ3Array)[i][j][0];

                j1FileRowMajor << j1Val << '\n';
                j2FileRowMajor << j2Val << '\n';
                j3FileRowMajor << j3Val << '\n';

            }
        }

        j1FileRowMajor.close();
        cout << "j1 saved" << "\n";
        j2FileRowMajor.close();
        cout << "j2 saved" << "\n";
        j3FileRowMajor.close();
        cout << "j3 saved" << "\n";
    }
}
/*==============================================================================
 * read from file - pareto
 *============================================================================*/
void paretoReadFromFile(multiarray *aJ1Vals, multiarray *aJ2Vals, multiarray *aJ3Vals, short int aNumPointsADim, short int aNumPointsBDim){
    /*
     Function: paretoReadFromFile

     Arguments: multiarrays of J1, J2, and J3 values at each (A,B) node on the Pareto grid, and
     number of points in the A and B dimensions.

     Purpose: Reads in values from the Pareto grid computation for use in finding
     the J3-Constrained Pareto front in Example 1.
     */

    ifstream j1FileRowMajor;
    ifstream j2FileRowMajor;
    ifstream j3FileRowMajor;

    if(gUncertainPareto){
        j1FileRowMajor.open("./output/j1_grid_uncert_row_maj.txt");
        j2FileRowMajor.open("./output/j2_grid_uncert_row_maj.txt");
        j3FileRowMajor.open("./output/j3_grid_uncert_row_maj.txt");
    }

    else if (gStationaryPareto){
        j1FileRowMajor.open("./output/j1_grid_stationary_row_maj.txt");
        j2FileRowMajor.open("./output/j2_grid_stationary_row_maj.txt");
        j3FileRowMajor.open("./output/j3_grid_stationary_row_maj.txt");
    }

    else if (gYRGPareto){
        j1FileRowMajor.open("./output/j1_grid_yrg_row_maj.txt");
        j2FileRowMajor.open("./output/j2_grid_yrg_row_maj.txt");
        j3FileRowMajor.open("./output/j3_grid_yrg_row_maj.txt");
    }

    short int timeEndIdx = 0;

    for(short int k = 0; k < timeEndIdx + 1; k++) {           // time
         for(short int i = 0; i < aNumPointsADim; i++) {      // distance
             for(short int j = 0; j < aNumPointsBDim; j++) {  // velocity
                 string j1Val;
                 string j2Val;
                 string j3Val;

                 getline(j1FileRowMajor, j1Val, '\n');
                 getline(j2FileRowMajor, j2Val, '\n');
                 getline(j3FileRowMajor, j3Val, '\n');

                 //convert string to double
                 double j1Double = stod(j1Val);
                 double j2Double = stod(j2Val);
                 double j3Double = stod(j3Val);

                 (*aJ1Vals)[i][j][k] = j1Double;
                 (*aJ2Vals)[i][j][k] = j2Double;
                 (*aJ3Vals)[i][j][k] = j3Double;

             } //end vel loop
         } //end pos loop
    }//end time loop
    j1FileRowMajor.close();
    cout << "j1 read" << "\n";
    j2FileRowMajor.close();
    cout << "j2 read" << "\n";
    j3FileRowMajor.close();
    cout << "j3 read" << "\n";

}


/*==============================================================================
 * Write to file - debug arrays - not in use
 *============================================================================*/
void debugWriteToFile (vector<double> aControlVec, vector<double> aVFVec){
    /*
     Function: trajWriteToFile

     Purpose: Writes results from trajectory tracer to file
     */

    /*---------------------------------Write to file -------------------------------------------------*/
    ofstream tfile;

    tfile.open("debug_negative_controls_info.txt"); //optimal trajectory
    //Optimal Trajectory write to file

    for(int i = 0;  i < gNumControls; i++) {
        //tfile << gDMax -  positionVals[pos];
        tfile <<  aControlVec[i];
      if(i != gNumControls - 1) {
        tfile << " ";
      }
    }
    tfile << endl;

    for(int j = 0; j < gNumControls; j++) {
      tfile << aVFVec[j];
        //cout << aVFVec[j] << "\n";
        if(j != gNumControls - 1) {
        tfile << " ";
      }
    }
    tfile << endl;
    for(int j = 0; j < gNumControls; j++) {
      tfile << gTerminalCondValsAtTestPoint[j];
        cout << gTerminalCondValsAtTestPoint[j] << "\n";
        if(j != gNumControls - 1) {
        tfile << " ";
      }
    }

    tfile << endl;
    for(int j = 0; j < gNumControls; j++) {
      tfile << gUNextValsAtTestPoint[j];
        //cout << velocityVals[vel] << "\n";
        if(j != gNumControls - 1) {
        tfile << " ";
      }
    }

    tfile << endl;
    for(int j = 0; j <= gNt - 1; j++) {
      tfile << gProbVals[j];
        //cout << velocityVals[vel] << "\n";
        if(j != gNt - 1) {
        tfile << " ";
      }
    }
    tfile << endl;
    for(int j = 0; j <= gNt - 1 ; j++) {
      tfile << (gNt - j) * gDT;
        //cout << velocityVals[vel] << "\n";
        if(j != gNt - 1) {
        tfile << " ";
      }
    }
    /*tfile << endl;
    for(int j = 0; j < gNumParabPts; j++) {
      tfile << gTestParabYPts[j];
        //cout << velocityVals[vel] << "\n";
        if(j != gNumParabPts - 1) {
        tfile << " ";
      }
    }

    tfile << endl;
    for(int j = 0; j < 4; j++) {
      tfile << gParabTestXNodes[j];
        //cout << velocityVals[vel] << "\n";
        if(j != 3) {
        tfile << " ";
      }
    }
    tfile << endl;

    for(int j = 0; j < 4; j++) {
      tfile << gParabTestYNodes[j];
        //cout << velocityVals[vel] << "\n";
        if(j != 3) {
        tfile << " ";
      }
    }*/

    tfile.close();
}
