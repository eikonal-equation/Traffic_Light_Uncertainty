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
 * File: InterpolationFunctions.cpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains all of the functions for imterpolation
 * routines used in the scheme to solve the HJB equations described in
 * the Optimal Driving manuscript and the control-interpolation used
 * while tracing the optimal trajectories.
 *============================================================================*/
#include "InterpolationFunctions.hpp"

//-----------------------------Libraries-----------------------------
#include <iostream>
#include <fstream>
#include <cmath>
#include "boost/multi_array.hpp"

//----------------------Project specific headers-----------------------
#include "variables.h"
#include "InitFunctions.hpp"
#include "CoordinateFunctions.hpp"
#include "ProbabilityFunctions.hpp"
#include "HelperFunctions.hpp"

/*==============================================================================
 * Bilinear Interpolation scheme
 *============================================================================*/
double bilinearInterp(multiarray *aU, multiarray * aOC, short int aTimeIndex, double aPos, double aVel, double aTime, double aControl){
    /*
     Function: bilinearInterp
     
     Arguments: value function array (aU), optimal control array, position (aPos), velocity (aVel),
     time (aTime), and control value
     
     Purpose: Carry out bilinear interpolation using four grid cell corners with finite value function value.
     */
    
    double interpValue = 0;
    
    //Convert physical coordinates to grid coordinates
    double dRescaled = (aPos - gDTarget) / gDPos;
    double vRescaled = aVel / gDVel;
    
    double tol = 1e-13;
    
    if (fabs(vRescaled) <= tol){
        vRescaled = 0;
    }
    
    //Cell corner grid coordinates
    short int dC = short(int(ceil(dRescaled)));
    short int dF = short(int(floor(dRescaled)));
    short int vC = short(int(ceil(vRescaled)));
    short int vF = short(int(floor(vRescaled)));
    
    //Determine corresponding physical coordinates
    double dCeiling = dC * gDPos + gDTarget;
    double dFloor = dF * gDPos + gDTarget;
    double vCeiling = vC * gDVel;
    double vFloor = vF * gDVel;
    
    //debugging
    if ((dF < 0)){
        dF = 0;
        dFloor = 0;
    }
    
    if (vF < 0){
        vF = 0;
        vFloor = 0;
    }
    
    if ((aPos > gDMax) || (aPos < gDTarget) || (aVel > gVMax) || (aVel < 0)){
        interpValue = gInfty;
    }
    
    else{
        assert(dC >= 0);
        assert(dF >= 0);
        assert(vF >= 0);
        assert(vC >= 0);
        assert(vC <= gVNum);
        assert(vF <= gVNum);
        assert(dRescaled >= 0);
        assert(vRescaled >= 0);
        
        assert(dRescaled >= 0);
        assert(vRescaled >= 0);
        
        
        //Value function values at each corner
        double uA = (*aU)[dF][vC][aTimeIndex];
        double uB = (*aU)[dC][vC][aTimeIndex];
        double uC = (*aU)[dF][vF][aTimeIndex];
        double uD = (*aU)[dC][vF][aTimeIndex];
        
        if (gControlInterp){
            uA = (*aOC)[dF][vC][aTimeIndex];
            uB = (*aOC)[dC][vC][aTimeIndex];
            uC = (*aOC)[dF][vF][aTimeIndex];
            uD = (*aOC)[dC][vF][aTimeIndex];
        }
        
        double beta = vRescaled - vF;
        double gamma = dRescaled - dF;
        
        
        assert(((1-beta)+beta)==1);
        assert(((1-gamma)+gamma)==1);
        
        //vertical interpolation
        
        double Q2 = (1-beta)*uC + beta*uA;
        double Q1 = (1-beta)*uD + beta*uB;
        
        interpValue = gamma*Q1 + (1-gamma)*Q2; //horizontal interpolated solution
    
        
        if(!gControlInterp){
            assert(interpValue >= 0);
        }
        
        
    }//end else
    return  interpValue;
    
}

/*======================================================================================
 * Cut cell bilinear Interpolation scheme for problems with max-accel bd
 *=====================================================================================*/
double maxAccelCutCellBilinearInterp (multiarray *aU, multiarray * aOC, short int aTimeIndex, double aNextPos, double aNextVel, double aNextTime, double aCurrentPos, double aCurrentVel, double aCurrentTime, double aControl){
    /*
     Function: maxAccelCutCellBilinearInterp()
     
     Arguments: value function array, optimal control array, index of timeslice to interpolate in,
     (d,v,t) at foot of characteristic (aNextPos, aNextVel, aNextTime), current (d,v,t), and control value
     
     Purpose: Bilinear interpolation to find value function at point when cell is cut by max-accel boundary.
     Utilizes cut-cells along d_{\beta} and appropriate boundary points to carry out the bilinear interpolation.
     
     * Interpolation cell:
     *
     * A-------B
     * |       |
     * Q1  p  Q2
     * |       |
     * C-------D
     *
     * Points:
     * A = (dFloor, vCeiling)
     * B = (dCeiling, vCeiling)
     * C = (dFloor, vFloor)
     * D = (dCeiling, vFloor)
     * p = (aD,aV)
     
     */
    
    /*---------------------------Grid Coordinate Init-------------------------------------------------*/
    double interpVal = 0;
    
    //Grid coordinates
    double dRescaled = (aNextPos - gDTarget) / gDPos;
    double vRescaled = aNextVel / gDVel;
    
    //Grid cell corners
    short int dC = short(int(ceil(dRescaled)));
    short int dF = short(int(floor(dRescaled)));
    short int vC = short(int(ceil(vRescaled)));
    short int vF = short(int(floor(vRescaled)));
    
    //Determine corresponding physical coordinates
    double dCeiling = dC * gDPos + gDTarget;
    double dFloor = dF * gDPos + gDTarget;
    double vCeiling = vC * gDVel;
    double vFloor = vF * gDVel;
    
    //debugging - catch if coords slightly out of bounds
    if ((dF < 0)){
        dF = 0;
        dFloor = 0;
    }
    
    if (vF < 0){
        vF = 0;
        vFloor = 0;
    }
    
    //Value function values at each corner
    double uA = (*aU)[dF][vC][aTimeIndex];
    double uB = (*aU)[dC][vC][aTimeIndex];
    double uC = (*aU)[dF][vF][aTimeIndex];
    double uD = (*aU)[dC][vF][aTimeIndex];
    
    if(gControlInterp){
        uA = (*aOC)[dF][vC][aTimeIndex];
        uB = (*aOC)[dC][vC][aTimeIndex];
        uC = (*aOC)[dF][vF][aTimeIndex];
        uD = (*aOC)[dC][vF][aTimeIndex];
    }
    
    //debugging
    bool isAllowedA = isAllowedRedPhase(vCeiling, dFloor, aNextTime);
    bool isAllowedB = isAllowedRedPhase(vCeiling, dCeiling, aNextTime);
    bool isAllowedC = isAllowedRedPhase(vFloor, dFloor, aNextTime);
    bool isAllowedP = isAllowed(aNextVel, aNextPos, aNextTime);
    bool isAllowedD = isAllowed(vFloor, dCeiling, aNextTime);
    
    double dParab = pow(aNextVel, 2) / (2* gAlpha);
    double timeToRedPhase = gTR - aNextTime;
    double timeToMaxV = (gVMax - aNextVel) / gBeta;
    double vCritical = gVMax - gBeta * timeToRedPhase; //velocity above which you'll reach vMax before TR if max accel
    double tol = gDPos / 10000;
    //Quad form for velocity at intersection between max braking and max accel parabs
    double parabExtensionA = (gBeta + gAlpha) / (2 * gAlpha * gBeta); //a
    double parabExtensionB = -gVMax / gBeta; //b
    double parabExtensionC = pow(gVMax, 2) / (2 * gBeta) -timeToRedPhase * gVMax; //c
    double parabExtensionVelThresh = (-parabExtensionB + sqrt(pow(parabExtensionB, 2) - 4 *parabExtensionA * parabExtensionC)) / (2 * parabExtensionA); //quadratic formula

    double velocityThreshRoot = gAlpha * timeToRedPhase * (1 + sqrt(1 + (gBeta / gAlpha))); //velocity at intersection between max braking parab and max accel line
    double velocityThresh = (vCritical < velocityThreshRoot) ? parabExtensionVelThresh : velocityThreshRoot; //set velocity threshold
    velocityThresh = (timeToRedPhase > 0) ? velocityThresh : -tol;
    
    
    /*---------------------------------Interpolation Scheme--------------------------------------------------------------*/
    
    //*******************************************Case One: one corner out***********************************************************
    if((isAllowedA) && (isAllowedB) && (isAllowedC) && (!isAllowedD)){
        double timeToRed = gTR - aNextTime;
        //Determine vBound
        bool nextPtIsSubjectToExtension = isSubjectToMaxAccelExtension(aNextVel, aNextPos, aNextTime);
        double vBound;
        double dBound;
        if (nextPtIsSubjectToExtension){
            double a = - 1 / (2 * gBeta);
            double b = gVMax / gBeta;
            double c = timeToRed * gVMax - (pow(gVMax, 2) / (2 * gBeta)) - dCeiling;
            vBound = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
            
            dBound =  ((-pow(aNextVel, 2)) / (2 * gBeta)) + ((aNextVel * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeToRed * gVMax);
            assert(dBound >= 0);
            assert(dBound >= dFloor);
            
            assert(dBound >= 0);
            assert(vBound >= 0);
        }
        
        else{
            vBound = (dCeiling / timeToRed) - 0.5 * gBeta * timeToRed;
            dBound = (aNextVel * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2));
            
            assert(dBound >= 0);
            assert(vBound >= 0);
    
        }
        
        double Q1;
        double Q2;
        double posWeight;
        double velWeight;
        
        //If vBound approx = vFloor, set equal
        double vFloorDiff = fabs(vBound - vFloor);
        double dCeilDiff = fabs(dBound - dCeiling);

        vBound = (vFloorDiff <= 1e-13) ? vFloor : vBound;
        assert(vBound >=0);
        assert(vBound >= vFloor);
        assert(vBound <= gVMax);
        
        assert(vBound <= vCeiling);
        
        //Interpolation between pw boundary and B not necessary to get Q1
        if (aNextVel <= vBound){
            /* Example cell:
             * Interpolation cell:
             *
             * A-------B
             * |      x| <-- (dBound, vBound)
             * Q2  p   |
             * | x     |
             * C-------D
             *
             * Points:
             * A = (dFloor, vCeiling) IN
             * B = (dCeiling, vCeiling) IN
             * C = (dFloor, vFloor) IN
             * D = (dCeiling, vFloor) OUT
             * p = (aD,aV) IN
             
             */
            
            assert(dBound >= dFloor);
            assert(dBound <= dCeiling);
            
            double dB = (dBound - gDTarget) / gDPos;
            Q1 = maxAccelBoundaryCost(aU, aOC, aNextVel, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
            
            if((gInTracer) && (gControlInterp)){
                Q1 = gBeta;
            }
            
            velWeight = vRescaled - vF;
            posWeight = (dRescaled - dF) / (dB - dF);
            assert(posWeight <= 1);
            assert(posWeight >= 0);
            assert(velWeight <= 1);
            assert(velWeight >= 0);
            
            assert(posWeight + (1 - posWeight) == 1);
            assert(velWeight + (1 - velWeight) == 1);
            
            Q2 = (1 - velWeight) * uC + velWeight * uA;
            
            interpVal = (1 - posWeight) * Q2 + posWeight * Q1;
        }
        
        //Interpolation necessary between boundary and B to get Q1
        else{
            //Vertical interp between B and bd to get Q1
            double vB = vBound / gDVel;
            double velWeight1 = (vRescaled - vB) / (vC - vB);
            //Catch case when vBound == vFloor
            if (vBound == vFloor){
                velWeight1 = vRescaled - vF;
            }
                        
            double rightBdPt = maxAccelBoundaryCost(aU, aOC, vBound, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
            
            if((gInTracer) && (gControlInterp)){
                rightBdPt = gBeta;
            }
            
            velWeight = vRescaled - vF;
            posWeight = dRescaled - dF;
            assert(posWeight <= 1);
            assert(posWeight >= 0);
            assert(velWeight <= 1);
            assert(velWeight >= 0);
            Q1 = (1 - velWeight1) * rightBdPt + velWeight1 * uB;
            
            Q2 = (1 - velWeight) * uC + velWeight * uA;
            
            interpVal = (1 - posWeight) * Q2 + posWeight * Q1;
            
            if(!gControlInterp){
                assert(interpVal >= 0);
            }
        }
        
    }//end one out
    
    //*******************************************Case two: two corners out****************************************************
    else if ((isAllowedA) && (isAllowedC) && (!isAllowedB) && (!isAllowedD)){
        bool nextPtIsSubjectToExtension = isSubjectToMaxAccelExtension(aNextVel, aNextPos, aNextTime);
        double Q1;
        double Q2;
        double posWeight;
        double velWeight;
        double dBound;
        
        if (nextPtIsSubjectToExtension){
            assert(timeToRedPhase > timeToMaxV);
            dBound =  ((-pow(aNextVel, 2)) / (2 * gBeta)) + ((aNextVel * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeToRedPhase * gVMax);
            assert(dBound >= 0);
        }
        
        else{
            dBound = (aNextVel * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2));
            assert(dBound >= 0);
        }
        
        assert(dBound >= dFloor);
        assert(dBound <= dCeiling);
        
        double dB = (dBound - gDTarget) / gDPos;
        posWeight = (dRescaled - dF) / (dB - dF);
        velWeight = vRescaled - vF;
        
        assert(posWeight <= 1);
        assert(posWeight >= 0);
        assert(velWeight <= 1);
        assert(velWeight >= 0);
        assert(posWeight + (1 - posWeight) == 1);
        assert(velWeight + (1 - velWeight) == 1);
    
        
        Q1 = maxAccelBoundaryCost(aU, aOC, aNextVel, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
        
        if((gInTracer) && (gControlInterp)){
            Q1 = gBeta;
        }
        
        Q2 = (1 - velWeight) * uC + velWeight * uA;
        
        interpVal = (1 - posWeight) * Q2 + posWeight * Q1;
        if(!gControlInterp){
            assert(interpVal >= 0);
        }
        
    }//end two out
        
    else if ((isAllowedA) && (isAllowedC) && (isAllowedD) && (!isAllowedB)){
            //Vert interp between A and C
            double beta = vRescaled - vF; //velocity weight
            double gamma = dRescaled - dF; //position weight
            assert(beta >= 0);
            assert(gamma >= 0);
            assert(((1-beta)+beta)==1);
            assert(((1-gamma)+gamma)==1);
            
            double Q2 = (1-beta)*uC + beta*uA;
            assert(Q2 >= 0);
            //Vert interp between parab and D
            double vBound;
            //Step 1: Determine Q1 and vBound
            double vBoundParab = sqrt(2 * dCeiling * gAlpha);
            vBound = vBoundParab;
            double vB = vBound / gDVel;
            double beta2 = (vRescaled - vF) / (vB - vF);
            
            //get d at bd intersection point - if dNext < threshold d, then use max accel pw, if not, use parabolic
            
            bool isSubjectToMaxAccelBD = isInAccelRegion(aNextVel, aNextPos, aNextTime);
            
            if(isSubjectToMaxAccelBD){
                //check which portion its subject to
                bool isSubjectToExtension = isSubjectToMaxAccelExtension(aNextVel, aNextPos, aNextTime);
                
                if(isSubjectToExtension){
                    assert(timeToRedPhase > timeToMaxV);
                    double dBound = ((-pow(aNextVel, 2)) / (2 * gBeta)) + ((aNextVel * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeToRedPhase * gVMax);
                    assert(dBound >= 0);
                    double dB = (dBound - gDTarget) / gDPos;
                    double gamma2 = (dRescaled - dF) / (dB - dF);
                    
                    assert(gamma2 <= 1);
                    assert(gamma2 >= 0);
                    
                    double Q1 = maxAccelBoundaryCost(aU, aOC, aNextVel, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
                    
                    if((gInTracer) && (gControlInterp)){
                        Q1 = gBeta;
                    }
                    
                    interpVal = gamma2 * Q1 + (1 - gamma2) * Q2; //interpolated solution
                    if(!gControlInterp){
                        assert(interpVal >= 0);
                    }
                    
                }
                
                else{
                    double dBound = (aNextVel * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2));
                    assert(dBound >= 0);
                    double dB = (dBound - gDTarget) / gDPos;
                    double gamma2 = (dRescaled - dF) / (dB - dF);
                    
                    assert(gamma2 <= 1);
                    assert(gamma2 >= 0);
                    double Q1 = maxAccelBoundaryCost(aU, aOC, aNextVel, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
                    
                    if((gInTracer) && (gControlInterp)){
                        Q1 = gBeta;
                    }
                    
                    interpVal = gamma2 * Q1 + (1 - gamma2) * Q2; //interpolated solution
                    if(!gControlInterp){
                        assert(interpVal >= 0);
                    }
                }
                
            }
        
        else{
            double Q1 = beta2 * parabBoundaryCost(aU, aOC, vBound, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget) + (1 - beta2) * uD;
            
            if((gInTracer) && (gControlInterp)){
                Q1 = beta2 * (-gAlpha) + (1 - beta2) * uD;
            }
            
            assert(Q1 >= 0);
            interpVal = gamma * Q1 + (1 - gamma) * Q2; //interpolated solution
            if(!gControlInterp){
                assert(interpVal >= 0);
            }
        }
        
    }
    
    //***********************************Case three: three corners out**************************************************
    else{
        double Q1;
        double Q2;
        double posWeight; //
        double velWeight;
        double dBound;
        double vBound;
        bool nextPtIsSubjectToExtension = isSubjectToMaxAccelExtension(aNextVel, aNextPos, aNextTime);
        
        if (nextPtIsSubjectToExtension){
            assert(timeToRedPhase > timeToMaxV);
            double a = - 1 / (2 * gBeta);
            double b = gVMax / gBeta;
            double c = timeToRedPhase * gVMax - (pow(gVMax, 2) / (2 * gBeta)) - dFloor;
            vBound = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
            dBound =  ((-pow(aNextVel, 2)) / (2 * gBeta)) + ((aNextVel * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeToRedPhase * gVMax);
            assert(dBound >= 0);
        }
        
        else{
            vBound = (dFloor / timeToRedPhase) - 0.5 * gBeta * timeToRedPhase;
            dBound = (aNextVel * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2));
            assert(dBound >= 0);
        }
        
        double leftBdPt = maxAccelBoundaryCost(aU, aOC, vBound, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
        
        if((gInTracer) && (gControlInterp)){
            leftBdPt = gBeta;
        }
        
        double vB = vBound / gDVel;
        double dB = (dBound - gDTarget) / gDPos;
        posWeight = (dRescaled - dF) / (dB - dF);
        velWeight = (vRescaled - vB) / (vC - vB);
        
        assert(dBound >= dFloor);
        assert(dBound <= dCeiling);
        assert(vBound <= gVMax);
        assert(vBound >= vFloor);
        assert(vBound <= vCeiling);
        assert(vBound >= 0);
        assert(posWeight <= 1);
        assert(posWeight >= 0);
        assert(velWeight <= 1);
        assert(velWeight >= 0);

        
        Q1 = maxAccelBoundaryCost(aU, aOC, aNextVel, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
        
        if((gInTracer) && (gControlInterp)){
            Q1 = gBeta;
        }
        
        Q2 = (1 - velWeight) * leftBdPt + velWeight * uA;
        
        interpVal = (1 - posWeight) * Q2 + posWeight * Q1;
        
        if(!gControlInterp){
            assert(interpVal >= 0);
        }
    }//end three out
    
    if(!gControlInterp){
        assert(interpVal >= 0);
    }
    return interpVal;
}

/*==============================================================================
 * Cut cell bilinear Interpolation scheme for parabolic portion
 *============================================================================*/
double parabCutCellBilinearInterp (multiarray *aU, multiarray *aOC, short int aTimeIndex, double aNextPos, double aNextVel, double aNextTime, double aCurrentPos, double aCurrentVel, double aCurrentTime, double aControl){
    
    /*
     Function: parabCutCellBilinearInterp
     
     Arguments: value function array, optimal control array, index of timeslice to interpolate in,
     (d,v,t) at foot of characteristic (aNextPos, aNextVel, aNextTime), current (d,v,t), and control value
     
     Purpose: Handles bilinear interpolation on parabolic portion of boundary.
     Uses cut cells along d_{\alpha} as appropriate to carry out the bilinear interpolation.
     
     * Interpolation cell:
     *
     * A-------B
     * |       |
     * Q1  p  Q2
     * |       |
     * C-------D
     *
     * Points:
     * A = (dFloor, vCeiling)
     * B = (dCeiling, vCeiling)
     * C = (dFloor, vFloor)
     * D = (dCeiling, vFloor)
     * p = (aD,aV)
     
     */
    
    /*-------------------Grid Cell Coordinates Init------------------------------------------*/
    
    double interpValue = 0; //return value

    double tNext = aTimeIndex * gDT; //physical time
    double tol = gDPos / 100000;
    
    //Grid coords
    double dRescaled = (aNextPos - gDTarget) / gDPos;
    double vRescaled = aNextVel / gDVel;
    
    //Grid cell corners
    short int dC = short(int(ceil(dRescaled)));
    short int dF = short(int(floor(dRescaled)));
    short int vC = short(int(ceil(vRescaled)));
    short int vF = short(int(floor(vRescaled)));
    
    //Determine corresponding physical coordinates
    double dCeiling = dC * gDPos + gDTarget;
    double dFloor = dF * gDPos + gDTarget;
    double vCeiling = vC * gDVel;
    double vFloor = vF * gDVel;
    
    bool isInParab = isInParabola(aNextVel, aNextPos, aNextTime);
    bool pointAInParab = isInParabola(vCeiling, dFloor, aNextTime);
    bool pointBInParab = isInParabola(vCeiling, dCeiling, aNextTime);
    bool pointCInParab = isInParabola(vFloor, dFloor, aNextTime);
    bool pointDInParab = isInParabola(vFloor, dCeiling, aNextTime);
    /*-----------------------------------Interp Scheme--------------------------------------------------------------------*/
    //check if dNext, vNext out of bounds
    if ((aNextPos < gDTarget) || (aNextVel > gVMax) || (aNextVel < 0) || (aNextPos > gDMax)){
        interpValue = gInfty;
    }
    
    else if(!isInParab){
        interpValue = gInfty;
    }
    
    else if (aNextPos == gDTarget){
        interpValue = 0;
    }
    
    else{
        //debugging
        if ((dF < 0)){
            dF = 0;
            dFloor = 0;
        }
        
        if (vF < 0){
            vF = 0;
            vFloor = 0;
        }
        
        //Value function values at each corner
        double uA = (*aU)[dF][vC][aTimeIndex];
        double uB = (*aU)[dC][vC][aTimeIndex];
        double uC = (*aU)[dF][vF][aTimeIndex];
        double uD = (*aU)[dC][vF][aTimeIndex];
        
        //ASIDE FOR CONTROL INTERP
        if(gControlInterp){
            uA = (*aOC)[dF][vC][aTimeIndex];
            uB = (*aOC)[dC][vC][aTimeIndex];
            uC = (*aOC)[dF][vF][aTimeIndex];
            uD = (*aOC)[dC][vF][aTimeIndex];
        }
        
        //************************************************Four corners in admissible region*****************************************
        if ((pointAInParab) && (pointBInParab) && (pointCInParab)){

            //use grid indices
            double beta = vRescaled - vF; //velocity weight
            double gamma = dRescaled - dF; //position weight

            assert(((1-beta)+beta)==1);
            assert(((1-gamma)+gamma)==1);
            
            double Q2 = (1-beta)*uC + beta*uA;
            double Q1 = (1-beta)*uD + beta*uB;

            interpValue = gamma*Q1 + (1-gamma)*Q2; //interpolated solution
            assert(!isnan(interpValue));
        }

        //No corners in admissible region
        else if (!pointDInParab){

            interpValue = gInfty;
        }
        
        //*****************************3 corners in*********************************************************
        else if ((!pointAInParab) && (pointBInParab) && (pointCInParab) && (uD < gInfty)){
            double beta = vRescaled - vF;
            
            assert(((1-beta)+beta)==1);
            
            double Q1 = beta * uB + (1-beta) * uD;
            
            //Step 2: Determine if vertical interpolation with bd is necessary
            double Q2;
            double dTilde; //either dFloor, or dBound depending on if interp between bd and C necessary
            double dBound;
            
            dBound = pow(aNextVel, 2) / (2 * gAlpha); //point on parabola
            
            assert(dBound <= dCeiling);
            
            //Interpolate between boundary and C to get Q2
            if (dBound < dFloor){
                double vBoundParab = sqrt(2 * dFloor * gAlpha);
                double vBoundOfficial = vBoundParab; //set once we know which one to take
                
                assert(vBoundOfficial <= vCeiling);
                assert(vBoundOfficial >= vFloor);
                
                double vB = vBoundOfficial / gDVel; //rescaled boundary v
                double beta2 = (vRescaled - vF) / (vB - vF); //velocity weight
                
                assert(beta2 >= 0);
                assert(beta2 <= 1);
                //Catch situations where vB and vF are equal
                if (vBoundOfficial == vFloor){
                    Q2 = parabBoundaryCost(aU, aOC, vBoundOfficial, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
                    
                    if((gInTracer) && (gControlInterp)){
                        Q2 = -gAlpha;
                    }
                    
                }
                
                else{
                    Q2 = beta2 * parabBoundaryCost(aU, aOC, vBoundOfficial, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget) + (1 - beta2) * uC;
                    
                    if((gInTracer) && (gControlInterp)){
                        Q2 = beta2 * (-gAlpha) + (1 - beta2) * uC;
                    }
                    
                }
                
                dTilde = dFloor;
            }
            
            else{
                Q2 = parabBoundaryCost(aU, aOC, aNextVel, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
                if((gInTracer) && (gControlInterp)){
                    Q2 = -gAlpha;
                }
                dTilde = dBound;
            }
            
            //Step 3: Interpolate between Q2 and Q1
            double dT = (dTilde - gDTarget) / gDPos;

            double gamma = (dRescaled - dT) / (dC - dT);
            assert(gamma <= 1);
            assert(gamma >= 0);
            assert(dCeiling != dTilde);
            assert((gamma + (1-gamma)) == 1);
            
            interpValue = Q2 * (1-gamma) + gamma * Q1;
            
        }//end 3 corners in
        
        //****************************************************2 Corners in*********************************************************
        else if ((!pointAInParab) && (pointBInParab) && (!pointCInParab) && (uD < gInfty)){
            //Step 1: Interpolation between D and B
            double beta = vRescaled - vF;
            double Q1 = beta * uB + (1-beta) * uD;
            
            //Step 2: Interpolation between boundary and Q1
            double dBound;
            dBound = pow(aNextVel, 2) / (2 * gAlpha);
            double Q2 = parabBoundaryCost(aU, aOC, aNextVel, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
            
            if((gInTracer) && (gControlInterp)){
                Q2 = -gAlpha;
            }
            
            assert(dBound >= dFloor);
            assert(dBound <= dCeiling);
            double dB = (dBound - gDTarget) / gDPos;
            double gamma = (dRescaled - dB) / (dC - dB);
            assert(gamma <= 1);
            assert(gamma >= 0);
            interpValue = Q2 * (1-gamma) + gamma * Q1;
            
        } //end 2 corners in
        
        //**********************************Top two corners out******************************************************
        else if ((!pointAInParab) && (!pointBInParab) && (pointCInParab) && (pointDInParab)){
            double Q1;
            double vRightBound = sqrt(2 * gAlpha * dCeiling);
            double vRB = vRightBound / gDVel;
            
            double beta = (vRescaled - vF) / (vRB - vF);
            double rightBDPt = parabBoundaryCost(aU, aOC, vRightBound, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
            
            if((gInTracer) && (gControlInterp)){
                rightBDPt = -gAlpha;
            }
            
            Q1 = beta * rightBDPt + (1 - beta) * uD;
            
            //Step 2: Determine if vertical interpolation with bd is necessary
            double Q2;
            double dTilde; //either dFloor, or dBound depending on if interp between bd and C necessary
            double dBound;
            
            dBound = pow(aNextVel, 2) / (2 * gAlpha); //point on parabola
            
            assert(dBound <= dCeiling);
            
            //Interpolate between boundary and C to get Q2
            if (dBound < dFloor){
                double vBoundOfficial = sqrt(2 * dFloor * gAlpha);
                
                assert(vBoundOfficial <= vCeiling);
                assert(vBoundOfficial >= vFloor);
                
                double vB = vBoundOfficial / gDVel; //rescaled boundary v
                double beta2 = (vRescaled - vF) / (vB - vF); //velocity weight
                
                assert(beta2 >= 0);
                assert(beta2 <= 1);
                //Catch situations where vB and vF are equal
                if (vBoundOfficial == vFloor){
                    Q2 = parabBoundaryCost(aU, aOC, vBoundOfficial, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
                    if((gInTracer) && (gControlInterp)){
                        Q2 = -gAlpha;
                    }
                }
                
                else{
                    Q2 = beta2 * parabBoundaryCost(aU, aOC, vBoundOfficial, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget) + (1 - beta2) * uC;
                    
                    if((gInTracer) && (gControlInterp)){
                        Q2 = beta2 * (-gAlpha) + (1 - beta2) * uC;
                    }
                    
                }
                
                dTilde = dFloor;
            }
            
            else{
                Q2 = parabBoundaryCost(aU, aOC, aNextVel, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
                if((gInTracer) && (gControlInterp)){
                    Q2 = -gAlpha;
                }
                dTilde = dBound;
            }
            
            //if on gridline, no need to interpolate in d
            double dCeilDiff = fabs(aNextPos - dCeiling);
            if (dCeilDiff <= tol){
                interpValue = Q1;
            }
            
            //Step 3: Interpolate between Q2 and Q1
            else{
                double dT = (dTilde - gDTarget) / gDPos;

                double gamma = (dRescaled - dT) / (dC - dT);
                
                assert(gamma <= 1);
                assert(gamma >= 0);
                assert(dCeiling != dTilde);
                assert((gamma + (1-gamma)) == 1);
                
                interpValue = Q2 * (1-gamma) + gamma * Q1;
            }
        }
        
        //******************************3 corners out***********************************************************
        else{
            assert(!pointAInParab);
            if(pointBInParab){
                
            }
            assert(!pointBInParab);
            assert(!pointCInParab);

            double vBound;
            double Q1;
            
            double dBound = pow(aNextVel, 2) / (2 * gAlpha);
            double Q2;
            
            //Step 1: Determine Q1 and vBound
            double vBoundParab = sqrt(2 * dCeiling * gAlpha);
            vBound = vBoundParab;
            assert(vBound <= vCeiling);
            assert(vBound >= vFloor);
            assert(dBound >= dFloor);
            assert(dBound <= dCeiling);
            double vB = vBound / gDVel;
            double beta = (vRescaled - vF) / (vB - vF);
            
            Q1 = beta * parabBoundaryCost(aU, aOC, vBound, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget) + (1 - beta) * uD;
            
            if((gInTracer) && (gControlInterp)){
                Q1 = beta * (-gAlpha) + (1 - beta) * uD;
            }
            
            //Step 2: Determine Q2 on boundary
            
            Q2 = parabBoundaryCost(aU, aOC, aNextVel, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
            
            if((gInTracer) && (gControlInterp)){
                Q2 = -gAlpha;
            }
            
            //Step 3: Interpolate between Q1 and Q2 to get value function at next point
            double dB = (dBound - gDTarget) / gDPos;
            double gamma = (dRescaled - dB) / (dC - dB);
            assert(gamma <= 1);
            assert(gamma >= 0);
            assert((gamma + (1-gamma)) == 1); //off for debugging - rounding issues when aD on bd
            
            interpValue = Q1 * gamma + (1-gamma) * Q2;
            
        }//end 1 corner in
        
        assert(!isnan(interpValue));
    }
    
    if(!gControlInterp){
        assert(interpValue >= 0);
    }
    assert(!isnan(interpValue));
    return interpValue;
}

/*============================================================================================================================
* Discontinuity Line Interp
*==========================================================================================================================*/
double discontinuityLineInterp(multiarray *aU, multiarray *aOC, short int aTimeIndex, double aNextPos, double aNextVel, double aNextTime, double aCurrentPos, double aCurrentVel, double aCurrentTime, double aControl){
    /*
     Function: discontinuityLineInterp
     
     Arguments: value function array, optimal control array, index of timeslice to interpolate in,
     (d,v,t) at foot of characteristic (aNextPos, aNextVel, aNextTime), current (d,v,t), and control value
     
     Purpose: Carries out interpolation over the discontinuity line in the "legal" region.
     Two cases to consider:
     1) When the origin corner is cut by the discont line, but the next point is still to the left of the line.
     2) When the next point can't keep up with the discont line
     */
    
    double interpVal = 0;
    
    //Grid coordinates
    double dRescaled = (aNextPos - gDTarget) / gDPos;
    double vRescaled = aNextVel / gDVel;
    
    //Grid cell corners
    short int dC = short(int(ceil(dRescaled)));
    short int dF = short(int(floor(dRescaled)));
    short int vC = short(int(ceil(vRescaled)));
    short int vF = short(int(floor(vRescaled)));
    
    //Determine corresponding physical coordinates
    double dCeiling = dC * gDPos + gDTarget;
    double dFloor = dF * gDPos + gDTarget;
    double vCeiling = vC * gDVel;
    double vFloor = vF * gDVel;
    
    //debugging - catch if coords slightly out of bounds
    if ((dF < 0)){
        dF = 0;
        dFloor = 0;
    }
    
    if (vF < 0){
        vF = 0;
        vFloor = 0;
    }
    
    //Value function values at each corner
    double uA = (*aU)[dF][vC][aTimeIndex];
    double uB = (*aU)[dC][vC][aTimeIndex];
    double uC = (*aU)[dF][vF][aTimeIndex];
    double uD = (*aU)[dC][vF][aTimeIndex];
    
    if(gInTerminalPenalty){
        assert(aTimeIndex == gNt);
    }
    
    if(gControlInterp){
        uA = (*aOC)[dF][vC][aTimeIndex];
        uB = (*aOC)[dC][vC][aTimeIndex];
        uC = (*aOC)[dF][vF][aTimeIndex];
        uD = (*aOC)[dC][vF][aTimeIndex];
    }
    
    //Max accel line locations
    double maxAccelPosAtVNext = maxAccelLine(aNextVel, aNextTime);
    double timeToMaxV = (gVMax - aNextVel) / gBeta;
    double timeToRedPhase = gTR - aNextTime;
    //Check one: If dNext, vNext ends up to the right of the max accel line - should NOT be considering controls that don't allow us to keep up with the max accel line
    if (aNextPos > maxAccelPosAtVNext){
        interpVal = gInfty;
    }
    
    else{
        //Max accel line locations
        double maxAccelVFloor = maxAccelLine(vFloor, aNextTime);
        double maxAccelVCeil = maxAccelLine(vCeiling, aNextTime);
        double vCritical = gVMax - gBeta * (gTR - aNextTime); //critical velocity above which driver breaks speed limit before light turns red
        if(gInTerminalPenalty){
            assert((gTR - aNextTime) == gYellowDuration);
        }
        /*-------------------Interpolation Scheme--------------------------------------------------------------*/
        
        //***********************Case One: one corner out***********************************************************

        if((dCeiling < maxAccelVFloor)){
            interpVal = bilinearInterp(aU, aOC, aTimeIndex, aNextPos, aNextVel, aNextTime, aControl);
        }
        
        else if((dFloor < maxAccelVCeil) && (dCeiling < maxAccelVCeil) && (dFloor < maxAccelVFloor) && (dCeiling > maxAccelVFloor)){
            assert(dCeiling > maxAccelLine(vFloor, aNextTime));
            assert(dFloor < maxAccelVFloor);
            assert(dFloor < maxAccelVCeil);
            assert(dCeiling < maxAccelVCeil);
            //Determine vBound
            double vBound;
            double dBound;

            if (aNextVel > vCritical){
                double a = - 1 / (2 * gBeta);
                double b = gVMax / gBeta;
                double c = timeToRedPhase * gVMax - (pow(gVMax, 2) / (2 * gBeta)) - dCeiling;
                vBound = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
                dBound =  ((-pow(aNextVel, 2)) / (2 * gBeta)) + ((aNextVel * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeToRedPhase * gVMax);
                assert(dBound >= 0);
                assert(vBound >=0);
                assert(dBound >= dFloor);
            }
            
            else{
                vBound = (dCeiling / timeToRedPhase) - 0.5 * gBeta * timeToRedPhase;
                dBound = (aNextVel * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2));
                assert(dBound >= 0);
                assert(vBound >= 0);
            }
            
            double Q1;
            double Q2;
            double posWeight; //
            double velWeight;
            
            //If vBound approx = vFloor, set equal
            double vFloorDiff = fabs(vBound - vFloor);
            
            vBound = (vFloorDiff <= 1e-13) ? vFloor : vBound;
            
            //Interpolation between pw boundary and B not necessary to get Q1
            if (aNextVel <= vBound){
                /* Example cell:
                 * Interpolation cell:
                 *
                 * A-------B
                 * |      x| <-- (dBound, vBound)
                 * Q2  p   |
                 * | x     |
                 * C-------D
                 *
                 * Points:
                 * A = (dFloor, vCeiling) IN
                 * B = (dCeiling, vCeiling) IN
                 * C = (dFloor, vFloor) IN
                 * D = (dCeiling, vFloor) OUT
                 * p = (aD,aV) IN
                 
                 */
                
                double dB = (dBound - gDTarget) / gDPos;
                double Q1Parab = parabCutCellBilinearInterp(aU, aOC, aTimeIndex, dBound, aNextVel, aNextTime, aCurrentPos, aCurrentVel, aCurrentTime, aControl);
                double Q1MaxAccel = maxAccelBoundaryCost(aU, aOC, aNextVel, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
                Q1 = min(Q1Parab, Q1MaxAccel);
                
                if((gInTracer) && (gControlInterp)){
                    Q1 = (Q1Parab < Q1MaxAccel) ? -gAlpha : gBeta;
                }
                
                velWeight = vRescaled - vF;
                posWeight = (dRescaled - dF) / (dB - dF);
                
                Q2 = (1 - velWeight) * uC + velWeight * uA;
                
                interpVal = (1 - posWeight) * Q2 + posWeight * Q1;
            }
            
            //Interpolation necessary between boundary and B to get Q1
            else{
                //Vertical interp between B and bd to get Q1
                double velWeight1 = vRescaled - (vBound / gDVel);
                double vB = vBound / gDVel;
                velWeight1 = (vRescaled - vB) / (vC - vB);
                
                //Catch case when vBound == vFloor
                if (vBound == vFloor){
                    velWeight1 = vRescaled - vF;
                }
                
                double rightBdPtMaxAccel = maxAccelBoundaryCost(aU, aOC, vBound, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
                
                double rightBDPtParab = parabCutCellBilinearInterp(aU, aOC, aTimeIndex, dCeiling, vBound, aNextTime, aCurrentPos, aCurrentVel, aCurrentTime, aControl);
                double rightBdPt = min(rightBDPtParab, rightBdPtMaxAccel);
                
                if((gInTracer) && (gControlInterp)){
                    rightBdPt = (rightBDPtParab < rightBdPtMaxAccel) ? -gAlpha : gBeta;
                }
                
                velWeight = vRescaled - vF;
                posWeight = dRescaled - dF;
                assert(posWeight <= 1);
                assert(posWeight >= 0);
                assert(velWeight <= 1);
                assert(velWeight >= 0);
                Q1 = (1 - velWeight1) * rightBdPt + velWeight1 * uB;
                
                Q2 = (1 - velWeight) * uC + velWeight * uA;
                
                interpVal = (1 - posWeight) * Q2 + posWeight * Q1;

                if(!gControlInterp){
                    assert(interpVal >= 0);
                }
            }
            
        }//end one out

        //*********************************Case two: two corners out****************************************************
        else if((dFloor < maxAccelVCeil) && (dCeiling > maxAccelVCeil) && (dFloor < maxAccelVFloor) && (dCeiling > maxAccelVFloor)){
            assert(dFloor < maxAccelVFloor);
            assert(dFloor < maxAccelVCeil);
            assert(dCeiling > maxAccelVCeil);
            assert(dCeiling > maxAccelVFloor);
            double Q1;
            double Q2;
            double posWeight;
            double velWeight;
            double dBound;
                
            if (aNextVel > vCritical){
                assert(timeToRedPhase > timeToMaxV);
                dBound = ((-pow(aNextVel, 2)) / (2 * gBeta)) + ((aNextVel * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeToRedPhase * gVMax);
                assert(dBound >= 0);
            }
            
            else{
                dBound = (aNextVel * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2));
                assert(dBound >= 0);
            }
            
            double dB = (dBound - gDTarget) / gDPos;
            posWeight = (dRescaled - dF) / (dB - dF);
            velWeight = vRescaled - vF;
            
            assert(dBound >= dFloor);
            assert(dBound <= dCeiling);
            assert(posWeight <= 1);
            assert(posWeight >= 0);
            assert(velWeight <= 1);
            assert(velWeight >= 0);
            assert(posWeight + (1 - posWeight) == 1);
            assert(velWeight + (1 - velWeight) == 1);

            
            double Q1MaxAccel = maxAccelBoundaryCost(aU, aOC, aNextVel, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
            double Q1Parab = parabCutCellBilinearInterp(aU, aOC, aTimeIndex, dBound, aNextVel, aNextTime, aCurrentPos, aCurrentVel, aCurrentTime, aControl);
            Q1 = min(Q1MaxAccel, Q1Parab);
            
            if((gInTracer) && (gControlInterp)){
                Q1 = (Q1Parab < Q1MaxAccel) ? -gAlpha : gBeta;
            }
            
            Q2 = (1 - velWeight) * uC + velWeight * uA;
            
            interpVal = (1 - posWeight) * Q2 + posWeight * Q1;
            if(!gControlInterp){
                assert(interpVal >= 0);
            }
            
        }//end two out

        else{
            assert(dCeiling > maxAccelVCeil);
            assert(dCeiling > maxAccelVFloor);
            assert(dFloor < maxAccelVCeil);
            assert(dFloor > maxAccelVFloor);
            assert(dFloor < maxAccelVCeil);
            assert(dCeiling > maxAccelVCeil);
            assert(dCeiling > maxAccelVFloor);

            double Q1;
            double Q2;
            double posWeight; //
            double velWeight;
            double dBound;
            double vBound;
            
            if(aNextVel > vCritical){
                assert(timeToRedPhase > timeToMaxV);
                double a = - 1 / (2 * gBeta);
                double b = gVMax / gBeta;
                double c = timeToRedPhase * gVMax - (pow(gVMax, 2) / (2 * gBeta)) - dFloor;
                vBound = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
                dBound =  ((-pow(aNextVel, 2)) / (2 * gBeta)) + ((aNextVel * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeToRedPhase * gVMax);
                assert(dBound >= 0);
            }
            
            else{
                vBound = (dFloor / timeToRedPhase) - 0.5 * gBeta * timeToRedPhase;
                dBound = (aNextVel * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2));
                assert(dBound >= 0);
            }
            assert(dBound >= dFloor);
            assert(dBound <= dCeiling);
            assert(vBound >= 0);
            assert(vBound >=vFloor);
            assert(vBound <= vCeiling);
            assert(vBound <= gVMax);
            
            double leftBdPt = maxAccelBoundaryCost(aU, aOC, vBound, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
            double leftBdParab = parabCutCellBilinearInterp(aU, aOC, aTimeIndex, dFloor, vBound, aNextTime, aCurrentPos, aCurrentVel, aCurrentTime, aControl);
            double leftBdMaxAccel = maxAccelBoundaryCost(aU, aOC, vBound, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
            leftBdPt = min(leftBdParab, leftBdMaxAccel);
            
            if((gInTracer) && (gControlInterp)){
                leftBdPt = (leftBdParab < leftBdMaxAccel) ? -gAlpha : gBeta;
            }
            
            double vB = vBound / gDVel;
            double dB = (dBound - gDTarget) / gDPos;
            posWeight = (dRescaled - dF) / (dB - dF);
            velWeight = (vRescaled - vB) / (vC - vB);
            
            assert(posWeight <= 1);
            assert(posWeight >= 0);
            assert(velWeight <= 1);
            assert(velWeight >= 0);

            Q1 = maxAccelBoundaryCost(aU, aOC, aNextVel, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
            double Q1Parab = parabCutCellBilinearInterp(aU, aOC, aTimeIndex, dBound, aNextVel, aNextTime, aCurrentPos, aCurrentVel, aCurrentTime, aControl);
            double Q1MaxAccel = maxAccelBoundaryCost(aU, aOC, aNextVel, aNextTime, gTR, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
            Q1 = min(Q1Parab, Q1MaxAccel);
            
            if((gInTracer) && (gControlInterp)){
                Q1 = (Q1Parab < Q1MaxAccel) ? -gAlpha : gBeta;
            }
            
            Q2 = (1 - velWeight) * leftBdPt + velWeight * uA;
            
            interpVal = (1 - posWeight) * Q2 + posWeight * Q1;
            
            if(!gControlInterp){
                assert(interpVal >= 0);
            }
        }//end three out
        
    }//end else
    return interpVal;
}

/*==============================================================================
 * Entire interpolation routine
 *============================================================================*/
double generalInterp(multiarray *aU, multiarray * aOC, multiarray *aUParab, multiarray *aOCParab, short int aTimeIndex, double aNextPos, double aNextVel, double aNextTime, double aCurrentPos, double aCurrentVel, double aCurrentTime, double aControl){
    /*
     * Function: generalInterp
     *
     * Arguments: value function and optimal control arrays for first and second passes of
     * the correction method, index of timeslice to interpolate in,
     * (d,v,t) at foot of characteristic (aNextPos, aNextVel, aNextTime), current (d,v,t), and control value
     *
     * Purpose: Carries out the full interpolation scheme in distance and velocity to determine value
     * of the value function at point p. Makes a determination about the
     * corners of the cell in relation to the admissible region or max accel line (as appropriate) and then
     * implements the proper scheme (bilinear interp with four cell corners or cut cell interp)
     * to determine the value at p.
     *
     *
     * Interpolation cell:
     *
     * A-------B
     * |       |
     * Q1  p  Q2
     * |       |
     * C-------D
     *
     * Points:
     * A = (dFloor, vCeiling)
     * B = (dCeiling, vCeiling)
     * C = (dFloor, vFloor)
     * D = (dCeiling, vFloor)
     * p = (aD,aV)
     *
     */
    
    /*----------------------Grid Cell Coordinates Init--------------------------------------------------------------------*/
    
    double interpValue = 0; //return value

    double tNext = aNextTime; //this is confusing, perhaps consider changing the name / structure of t here
    
    //Convert physical coordinates to grid coordinates
    double dRescaled = (aNextPos - gDTarget) / gDPos;
    double vRescaled = aNextVel / gDVel;

    //Cell corner grid coordinates
    short int dC = short(int(ceil(dRescaled)));
    short int dF = short(int(floor(dRescaled)));
    short int vC = short(int(ceil(vRescaled)));
    short int vF = short(int(floor(vRescaled)));
    
    //Determine corresponding physical coordinates
    double dCeiling = dC * gDPos + gDTarget;
    double dFloor = dF * gDPos + gDTarget;
    double vCeiling = vC * gDVel;
    double vFloor = vF * gDVel;

/*--------------------------------------Interpolation Scheme--------------------------------------------------------------*/
    //check if dNext, vNext out of bounds
    if ((aNextPos < gDTarget) || (aNextVel > gVMax) || (aNextVel < 0) || (aNextPos > gDMax)){
        interpValue = gInfty;
    }
    
    else if (aNextPos == gDTarget){
        interpValue = 0;
    }
    
    else{
        //debugging
        if ((dF < 0)){
            dF = 0;
            dFloor = 0;
        }
        
        if (vF < 0){
            vF = 0;
            vFloor = 0;
            assert(aCurrentVel == 0);
        }

        double maxAccelLinePos = maxAccelLine(aNextVel, aNextTime);
        double maxAccelLineCurrentPos = maxAccelLine(aCurrentVel, aCurrentTime);
        
        double dLeft = dFloor;
        double vUpper = vCeiling;
        
        bool isAllowedTopLeftStencilCorner = isAllowed(vUpper, dLeft, aNextTime);
        if((gGreenRedGreenProb) && (gFirstPass)){
            bool isInParab = isInParabola(aNextVel, aNextPos, aNextTime);
            double dParab = pow(aNextVel, 2) / (2 * gAlpha);

            if((!isInParab) && (aNextPos > 0) && (aNextTime < gTG)){
                interpValue = gInfty;
            }

            else if ((aNextPos > 0) && (aNextPos == dParab) && (aNextTime < gTG)){
                interpValue = parabBoundaryCost(aU, aOC, aNextVel, aNextTime, gTG, gGamma, gC1, gC2, gC3, gDPos, gDVel, gDTarget);
                assert(!isnan(interpValue));
                assert(interpValue >= 0);
            }

            else if((aNextPos > 0) && (aNextTime < gTG)){
                interpValue = parabCutCellBilinearInterp(aU, aOC, aTimeIndex, aNextPos, aNextVel, aNextTime, aCurrentPos, aCurrentVel, aCurrentTime, aControl);

                assert(!isnan(interpValue));
                if(!gInTracer){
                    assert(interpValue >= 0);
                }
            }

            else{
                interpValue = bilinearInterp(aU, aOC, aTimeIndex, aNextPos, aNextVel, aNextTime, aControl);
                assert(!isnan(interpValue));
                if(!gInTracer){
                    assert(interpValue >= 0);
                }
            }
        }

        else if ((gGreenRedGreenProb) && (!gFirstPass)){
            if(((aCurrentPos <= maxAccelLineCurrentPos) && (aNextPos > maxAccelLinePos))){
                interpValue = gInfty;
            }
            
            else if ((gInTerminalPenalty) && (aNextPos > maxAccelLinePos)){
                //BL Interp with first pass YRG
                interpValue = bilinearInterp(aUParab, aOCParab, 0, aNextPos, aNextVel, gTerminalT, aControl);
            }
            
            else{
                
                if ((aNextPos > 0) && (aNextPos <= maxAccelLinePos)){
                    interpValue = discontinuityLineInterp(aU, aOC, aTimeIndex, aNextPos, aNextVel, aNextTime, aCurrentPos, aCurrentVel, aCurrentTime, aControl);
                }//end max accel line interp
                
                else{
                    double dParab = pow(aNextVel, 2) / (2 * gAlpha);
                    
                    if(aNextPos >= dParab){
                        interpValue = parabCutCellBilinearInterp(aU, aOC, aTimeIndex, aNextPos, aNextVel, aNextTime, aCurrentPos, aCurrentVel, aCurrentTime, aControl);
                    }
                    else{
                        interpValue = bilinearInterp(aU, aOC, aTimeIndex, aNextPos, aNextVel, aNextTime, aControl);
                    }
                }//end not max accel else
            }//end legal else
        }//end grg first pass else

        else{
            if(gInTerminalPenalty){
                double currentT = gTerminalT - gDT;
                double parabOnlyInterpValue = parabCutCellBilinearInterp(aUParab, aOCParab, 0, aNextPos, aNextVel, gTerminalT, aCurrentPos, aCurrentVel, currentT, aControl);
                double discontLineInterpVal = 0;
                if(((aCurrentPos <= maxAccelLineCurrentPos) && (aNextPos > maxAccelLinePos))){
                    interpValue = gInfty;
                }
                
                else if (aNextPos > maxAccelLinePos){
                    //BL Interp with first pass YRG
                    interpValue = bilinearInterp(aUParab, aOCParab, 0, aNextPos, aNextVel, gTerminalT, aControl);
                }
                else if ((aNextPos > 0) && (aNextPos <= maxAccelLinePos) && (aCurrentPos <= maxAccelLineCurrentPos)){
                    discontLineInterpVal = discontinuityLineInterp(aU, aOC, aTimeIndex, aNextPos, aNextVel, aNextTime, aCurrentPos, aCurrentVel, aCurrentTime, aControl);
                    interpValue = min(parabOnlyInterpValue, discontLineInterpVal);
                }
                else{
                    interpValue = bilinearInterp(aU, aOC, aTimeIndex, aNextPos, aNextVel, aNextTime, aControl);
                }
                
            }
            
            else{
                interpValue = bilinearInterp(aU, aOC, aTimeIndex, aNextPos, aNextVel, aNextTime, aControl);
                assert(!isnan(interpValue));
            }
            if(!gControlInterp){
                assert(interpValue >= 0);
            }
        }
        
    }
    return interpValue;
}

/*==============================================================================
* Interpolation in d,v,t
*============================================================================*/
double dvtInterpolation(multiarray *aU, multiarray *aOC, multiarray *aUParab, multiarray *aOCParab, double aV, double aD, double aDistCurrent, double aVelCurrent, double aCurrentTime, double aTimeUpdated, double aControlValue){
    
    /*
     Function: dvtInterpolation
     
     Arguments: value function and optimal control arrays for first and second pass of
     interpolation scheme, (d,v) at foot of characteristic (aD, aV),
     current (d,v,t), time at foot of characteristic (aTimeUpdated) and control value
     
     Purpose: Carries out full trilinear interpolation in d,v,t. Interpolates bilinearly
     in upper and lower timeslices. Then interpolates linearly between those two results
     to obtain value.
     
     */
    
    double interpVal;
    
    short int lowerTimeIdx = floor(aTimeUpdated / gDT);
    short int upperTimeIdx = ceil(aTimeUpdated / gDT);
    
    assert(lowerTimeIdx >= 0);
    assert(upperTimeIdx <= gNt);
    
    double lowerTime = lowerTimeIdx * gDT;
    double upperTime = upperTimeIdx * gDT;
    
    assert(lowerTime >= 0);
    assert(upperTime <= gTerminalT);
    
    double lowerInterp = generalInterp(aU, aOC, aUParab, aOCParab, lowerTimeIdx, aD, aV, lowerTime, aDistCurrent, aVelCurrent, aCurrentTime, aControlValue);
    double upperInterp = generalInterp(aU, aOC, aUParab, aOCParab, upperTimeIdx, aD, aV, upperTime, aDistCurrent, aVelCurrent, aCurrentTime, aControlValue);
    
    double timeWeight = (aTimeUpdated - lowerTime) / gDT;
    
    interpVal = timeWeight * upperInterp + (1 - timeWeight) * lowerInterp;
    
    if(!gControlInterp){
        assert(interpVal >= 0);
    }
    
    return interpVal;
}

/*==============================================================================
 * Interpolating between control values
 *============================================================================*/
struct sOptimalValuesControlInterp controlInterp(multiarray *aU, multiarray *aOC, multiarray *aUFirstPass, multiarray *aOCFirstPass, double aCurrentPos, double aCurrentVel, double aCurrentTime, short int aCurrentTimeIndex, double aTau, double aProbOfLightChange){
    /*
     Function: controlInterp
     
     Arguments: value function and optimal control arrays for first and second passes of
     the correction method, index of timeslice to interpolate in,
     (d,v,t) at foot of characteristic (aNextPos, aNextVel, aNextTime), current (d,v,t), and control value
     
     Purpose: Interpolates the optimal control values at the gridpoints determined in the PDE solve
     to find optimal control value in trajectory tracing. Uses the full bilinear interpolation scheme executed by
     generalInterp. After determining the optimal control at aCurrentTime, computes dNext, vNext using the control.
     
     Returns: struct containing optimal control value at (aCurrentPos, aCurrentVel, aCurrentTime) and
     (dNext, vNext, tNext) using this control.
     */
    sOptimalValuesControlInterp controlInterpVals;
    
    //Step 1: Determine the grid cell that (aCurrentPos, aCurrentVel) falls in
    double interpolatedControlValue = generalInterp(aU, aOC, aUFirstPass, aOCFirstPass, aCurrentTimeIndex, aCurrentPos, aCurrentVel, aCurrentTime, aCurrentPos, aCurrentVel, aCurrentTime, 0);
    
    double dNextOfficial = aCurrentPos - (aCurrentVel * aTau + 0.5 * interpolatedControlValue * pow(aTau, 2));
    double vNextOfficial = aCurrentVel + interpolatedControlValue * aTau;
    double tNextOfficial = aCurrentTime + aTau;
    
    controlInterpVals.dNextOfficial = dNextOfficial;
    controlInterpVals.vNextOfficial = vNextOfficial;
    controlInterpVals.tNextOfficial = tNextOfficial;
    controlInterpVals.oc = interpolatedControlValue;
    
    
    return controlInterpVals;
    
}

