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
 * File: CoordinateFunctions.cpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains the tau refinement function which adjusts SL step tau
 * when the vehicle will cross the boundary into the illegal region in an amount of time
 * less than \Delta t.
 *
 * This file also contains functions responsible for determining whether
 * a (d,v) point belongs to the allowable region at time t, where a (d,v) point
 * lies in relation to the pw boundary at time t, the location of the maximum
 * acceleration line d_{\beta} at time t, and the threshold velocity at time t.
 *
 *============================================================================*/

#include "CoordinateFunctions.hpp"

//-----------------------------Libraries-------------------------------------
#include <iostream>
#include <fstream>
#include <cmath>
#include "boost/multi_array.hpp"

//--------------------Project specific header files----------------------------
#include "variables.h"
#include "StationarySolver.hpp"
#include "ProblemFunctions.hpp"
#include "InitFunctions.hpp"

/*==============================================================================
* Determine dNext, vNext, tNext and tau refined
*============================================================================*/
struct sNextCoordinates tauRefinement(double aD, double aV, double aT, double aControl, double aTau){
    /*
     * Function: tauRefinement
     *
     * Arguments: position (aD), velocity (aV), time (aV), control value (aT),
     * and SL timestep (aTau)
     *
     * Purpose: Takes in current distance and velocity values and determines
     * the amount of time necessary to reach the boundary with current speed from
     * current position. If this time is less than the timestep size DT, then
     * there is no need to interpolate to find the value of the value function at
     * dNext, vNext, t+DT.
     *
     * Assumption: tauRefinement is run after setting value function of points that
     * start outside of legal region to infinity
     *
     * Relies heavily on isAllowed function, threshVelocityMaxAccel function, and
     * isOnAllowedPiece function to determine the point's location, and where the
     * point ends up relative to the pw boundary.
     *
     * Returns: A struct containing the value of tau used in the interpolation and dNext,
     * vNext, and tNext computed with this tau and the supplied control value (aControl)
     */
    
    sNextCoordinates sCoords;

    double tau = aTau;
    double epsilon = 1e-13;
    
    //compute next point
    double dNext = aD - (aV * aTau + 0.5 * aControl * pow(aTau, 2));
    double vNext = aV + aControl * aTau;
    double tNext = aT + aTau;
    
    double dParab = pow(aV, 2) / (2 * gAlpha); //parabolic braking boundary
    double threshV = threshVelocityMaxAccel(aD, aV, aT);
    
    //checking the location of the set of points above to determine if they are legal or illegal
    bool nextLocIsAdmissible = isAllowed(vNext, dNext, tNext); //check to see where point ends up in next timeslice
    bool startingLocIsOK = isAllowedRedPhase(aV, aD, aT);
    bool startsInAccelRegion = isInAccelRegion(aV, aD, aT);
    
    //Step 1: Determine if next location admissible or not. If not admissible, then point crossed bd in time less than tau.
    if ((gFirstPass) && (gGreenRedGreenProb)){
        double a = (pow(aControl, 2) / (2 * gAlpha)) + 0.5 * aControl;
        double b = ((aV * aControl) / gAlpha) + aV;
        double c = (pow(aV, 2) / (2 * gAlpha)) - aD;
        
        double discriminant = pow(b, 2) - (4 * a * c);
        double newTau = (-b + sqrt(discriminant)) / (2 * a);
        
        if ((aControl == 0) && (!nextLocIsAdmissible)){
            newTau = (aD / aV) - (aV / (2*gAlpha));
            tau = (newTau < aTau) ? newTau : aTau;
            double newD = aD - (aV * tau + 0.5 * aControl * pow(tau, 2));
            double newV = aV + aControl * tau;
            double newT = aT + tau;
            bool updatedPtOnParab = isOnParabPiece(newD, newV, newT);
            if(tau < gDT){
                assert(isAllowed(newV, newD, newT));
                assert(updatedPtOnParab);
            }
        }
        
        //set tau based on time to BD
        if ((aControl != -gAlpha) && (!isnan(newTau)) && (aT + newTau < gTG) && (!nextLocIsAdmissible) && (vNext <= gVMax)){
            tau = (newTau < aTau) ? newTau : aTau;
            assert(aTau <= gDT);
            assert(newTau >= 0);
            double newD = aD - (aV * tau + 0.5 * aControl * pow(tau, 2));
            double newV = aV + aControl * tau;
            double newT = aT + tau;
            bool updatedPtOnParab = isOnParabPiece(newD, newV, newT);
            bool updatedPtAllowed = isAllowed(newV, newD, newT);
            
            if(tau < gDT){
                if(!updatedPtOnParab){
                    bool updatedPtOnParabTest = isOnParabPiece(newD, newV, newT);
                }
                
                if(!updatedPtAllowed){
                    bool updatedPtAllowedTest = isAllowed(newV, newD, newT);
                }
                
                assert(updatedPtAllowed);
                assert(updatedPtOnParab);
            }
        }
        
        else{
            tau = aTau;
        }
    }
    
    else if ((!gFirstPass) && (aD > 0) && (aD > dParab) && (!nextLocIsAdmissible) && (aV >= threshV) && (gHasBoundaryOn)){
        double a = (pow(aControl, 2) / (2 * gAlpha)) + 0.5 * aControl;
        double b = ((aV * aControl) / gAlpha) + aV;
        double c = (pow(aV, 2) / (2 * gAlpha)) - aD;
        
        double discriminant = b * b - 4 * a * c;
        double newTau = (-b + sqrt(discriminant)) / (2 * a);
        
        if (aControl == 0){
            newTau = (aD / aV) - (aV / (2*gAlpha));
            tau = (newTau < aTau) ? newTau : aTau;
        }
        
        //set tau based on time to BD
        else if ((aControl != -gAlpha) && (!isnan(newTau)) && (aT + tau < gTG)){
            tau = (newTau < aTau) ? newTau : aTau;
        }
        
        else{
            tau = aTau;
        }
    }
    
    else{
        //Control value not >= 0 and/or time to reach boundary or target from starting position is >= gDT
        //Check to see if vehicle hits target in time less than DT
        if (dNext < gDTarget){
            assert(dNext <0);
            if (aControl != 0){
                double discriminant = pow(aV, 2) + 2 * aControl * (-gDTarget + aD);
            
                tau = (-aV + sqrt(discriminant)) / aControl;
            
                assert(tau <= aTau);
            }
            
            else if ((aControl == 0) && (aV != 0)){
                tau = (aD - gDTarget) / aV;
                assert(tau <= aTau);
            }
        }
        
    }
    
    //Small rounding correction for very small tau
    if (abs(tau) <= 10e-13){
        tau = 0;
    }
    
    if(gUncertainGreenProb){
        assert(tau == gDT);
    }
    
    //-------------------Compute and save dNext, vNext, tNext--------------------------
    double dNextOfficial = aD - (aV * tau + 0.5 * aControl * pow(tau, 2));
    double vNextOfficial = aV + aControl * tau;
    double tNextOfficial = aT + tau;
    
    dNextOfficial = ((tau < aTau) && (aD > 0)) ? pow(vNextOfficial, 2) / (2 * gAlpha) : dNextOfficial;
    dNextOfficial = ((tau < aTau) && (aD < 0)) ? gDTarget : dNextOfficial;
    
    sCoords.tauValue = tau;
    sCoords.nextPos = dNextOfficial;
    sCoords.nextVel = vNextOfficial;
    sCoords.nextTime = tNextOfficial;
    
    assert(tau >= 0);
    
    return sCoords;
}

/*==============================================================================
 * Determining if (d,v) point in allowable region at time t
 *============================================================================*/
bool isAllowed(double aV, double aD, double aT){
    /*
     * Function: isAllowed
     *
     * Arguments: velocity (aV), position (aD), time (aT)
     *
     * Purpose: determines the location of a specified point (aD,aV) in relation
     * to the boundary at time aT. The domain is partitioned into an admissible region and
     * an illegal region separated by a piecewise boundary that is only in effect
     * during red light and yellow light phases.
     *
     * Relies on the isAllowedRedPhase function to make the allowed/disallowed
     * determination during the yellow and red phases.
     *
     * Returns: True if the point belongs to the allowable
     * region at time t and false otherwise.
     */
    
    bool    allowed         = true;
    
    //--------------------break down into cases----------------------------------
    
    //Case 1: Out of bounds
    if ((aV < 0) || (aV > gVMax) || (aD > gDMax)) {
        allowed = false;
    }
    
    //Case 2: In bounds but in yellow or red phase, piecewise boundary in effect
    else if ((gGreenRedGreenProb) && (aT < gTG)){
        allowed = isAllowedRedPhase(aV, aD, aT);
    }
    
    return allowed;
}

/*==============================================================================
 * Determining if d,v in allowable region during the yellow and red phases
 *============================================================================*/
bool isAllowedRedPhase(double aV, double aD, double aT){
    /*
     * Function: isAllowedRedPhase
     *
     * Arguments: velocity (aV), position (aD), time (aT)
     *
     * Purpose: determines the location of a specified point (aD,aV) in relation
     * to the boundary at time aT. The domain is partitioned into an admissible region and
     * an illegal region separated by a piecewise boundary that is only in effect
     * during red light and yellow light phases.
     *
     * Returns: True if the point belongs to the allowable
     * region at time t during the yellow or red phase and false otherwise.
     */
    
    bool isAllowedRP = true;
    
    double timeToRedStart = gTR - aT; //time to red light
    double epsilon = 1e-7;
    double timeToRedPhase = (timeToRedStart > 0) ? timeToRedStart : 0;
    double timeToMaxV = (gVMax - aV) / gBeta; //time to reach max velocity
    
    double parabolicBound = (pow(aV,2) / (2 * gAlpha)); //max braking parab
    double maxAccelLine = (aV * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2)); //max accel line
    double maxAccelParabExtension = ((-pow(aV, 2)) / (2 * gBeta)) + ((aV * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeToRedPhase * gVMax); //parabolic extension when timeToMaxV < timeToRedPhase
    
    double vCritical = gVMax - gBeta * timeToRedPhase; //velocity above which you reach vMax before TR when accel maximally
    
    //quad formula for velocity at intersection between parab max accel extension and max braking parab
    double parabExtensionA = (gBeta + gAlpha) / (2 * gAlpha * gBeta); //a
    double parabExtensionB = -gVMax / gBeta; //b
    double parabExtensionC = pow(gVMax, 2) / (2 * gBeta) -timeToRedPhase * gVMax; //c
    double parabExtensionVelThresh = (-parabExtensionB + sqrt(pow(parabExtensionB, 2) - 4 *parabExtensionA * parabExtensionC)) / (2 * parabExtensionA); //quad form
    
    double velocityThreshRoot = gAlpha * timeToRedPhase * (1 + sqrt(1 + (gBeta / gAlpha))); //velocity at intersection between max accel line and max braking parab
    double velocityThresh = (vCritical < velocityThreshRoot) ? parabExtensionVelThresh : velocityThreshRoot; //set threshold v accordingly
    velocityThresh = (timeToRedPhase >= 0) ? velocityThresh : -epsilon;
    
    //First pass modification - no max accel line at first pass
    velocityThresh = (gFirstPass) ? -epsilon : velocityThresh;
    
    assert(velocityThresh >= -epsilon);
    
    //***********Determining if point is in illegal space between max accel region and braking parab*******
    
    //Out of bounds, valid in green, red, and uncertain green cases
    if ( (aV < 0) || (aV > gVMax) || (aD > gDMax)) {
        isAllowedRP = false;
    }
    
    else if ((velocityThresh <= 0) && (!gInUncertainGreenPhase) && (gFirstPass) && (aT < gTG)){
        if ((aD > 0) && (fabs(aD - parabolicBound) > epsilon) && (aD < parabolicBound)){
            isAllowedRP = false;
        }
    }
    
    //vCritical less than max accel v intersection with parabolic bd, only the parabolic accel extension
    else if (velocityThresh == parabExtensionVelThresh){
        if ((aD < parabolicBound) && (aV >= velocityThresh) && (aD > maxAccelParabExtension)){
            isAllowedRP = false;
        }
    }
    
    else if (velocityThresh == velocityThreshRoot){
        if ((aD > maxAccelLine) && (aV >= velocityThresh) && (aV <= vCritical) && (aD < parabolicBound)){
            isAllowedRP = false;
        }
        
        //Between parabolic extension to max accel line and parabolic bd
        else if ((aD < parabolicBound) && (aD > maxAccelParabExtension) && (aV > velocityThresh) && (aV >= vCritical)){
            isAllowedRP = false;
        }
    }
    
    else if (velocityThresh <= 0){
        if ((aD >= 0) && (aD < parabolicBound)){
            isAllowedRP = false;
        }
    }
    
    return isAllowedRP;
}

/*=====================================================================================
 * Determining if starting d,v is in the legal region above max accel line
 *=======================================================================================*/
bool isInAccelRegion(double aV, double aD, double aT){
    /*
     Function: isInAccelRegion
     
     Arguments: velocity (aV), position (aD), time (aT)
     
     Purpose: Determines if (d,v) point falls in the region to the LEFT of the maximum
     acceleration line d_{\beta}.
     
     Returns: True if the (d,v) point falls to the left of d_{\beta} and false otherwise.
     */
    
    bool inAccelRegion = false;
    
    if (gHasBoundaryOn){
        
        double timeToRedStart = gTR - aT; //time to red phase
        double epsilon = 1e-13;
        double timeToRedPhase = (timeToRedStart > 0) ? timeToRedStart : 0;
        double timeToMaxV = (gVMax - aV) / gBeta; //time to reach max velocity
        
        double parabolicBound = (pow(aV,2) / (2 * gAlpha)); //max braking parab
        double maxAccelLine = (aV * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2)); //max accel line
        double maxAccelParabExtension = ((-pow(aV, 2)) / (2 * gBeta)) + ((aV * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeToRedPhase * gVMax); //parabolic extension when timeToMaxV < timeToRedPhase
        double vCritical = gVMax - gBeta * timeToRedPhase; //velocity above which you'll reach vMax before TR if max accel
        
        //Quad form for velocity at intersection between max braking and max accel parabs
        double parabExtensionA = (gBeta + gAlpha) / (2 * gAlpha * gBeta); //a
        double parabExtensionB = -gVMax / gBeta; //b
        double parabExtensionC = pow(gVMax, 2) / (2 * gBeta) -timeToRedPhase * gVMax; //c
        double parabExtensionVelThresh = (-parabExtensionB + sqrt(pow(parabExtensionB, 2) - 4 *parabExtensionA * parabExtensionC)) / (2 * parabExtensionA); //quadratic formula

        double velocityThreshRoot = gAlpha * timeToRedPhase * (1 + sqrt(1 + (gBeta / gAlpha))); //velocity at intersection between max braking parab and max accel line
        double velocityThresh = (vCritical < velocityThreshRoot) ? parabExtensionVelThresh : velocityThreshRoot; //set velocity threshold
        velocityThresh = (timeToRedPhase > 0) ? velocityThresh : -epsilon;
        
        //***********************************************Determine if in the acceleration region above max braking parab************************
        
        //Between d= 0 and parabolic portion of the bd
        
        if (aT >= gTR){
            inAccelRegion = false;
        }
        
        else if ((aD >= 0) && (aD < parabolicBound) && (aV <= velocityThresh)){
            inAccelRegion = true;
        }
        
        //vCritical less than max accel line v intersection with parabolic bd, only the parabolic accel extension
        else if (velocityThresh == parabExtensionVelThresh){
            if ((aD >= 0) && (aV >= velocityThresh) && (aD <= maxAccelParabExtension)){
                inAccelRegion = true;
            }
        }
        
        //vCritical greater than max accel line velocity at intersection with max braking bd
        else if (velocityThresh == velocityThreshRoot){
            //Subject to max accel line
            if ((aD >= 0) && (aD <= maxAccelLine) && (aV >= velocityThresh) && (aV < vCritical)){
                inAccelRegion = true;
            }
            
            //Subject to parabolic-max-accel extension
            else if ((aD >= 0) && (aD <= maxAccelParabExtension) && (aV > velocityThresh) && (aV >= vCritical)){
                inAccelRegion = true;
            }
        }
        
        return inAccelRegion;
    }//end in red phase
    return inAccelRegion;
}

/*=======================================================================================
 * Location of (d,v) point in relation to PARABOLIC boundary d_{\alpha}
 *=====================================================================================*/
bool isInParabola(double aV, double aD, double aT){
    /*
     Function: isInParabola
     
     Arguments: velocity (aV), position (aD), time (aT)
     
     Purpose: Determines which side of parabolic boundary vehicle is on / if possible to come to a complete stop before light change with max braking
     
     Returns: True if ON or IN parab only!
     */
    
    bool isAllowedParab = false;
    
    double parabolicBound = (pow(aV,2) / (2 * gAlpha)); //parabolic braking bd
    
    if ( (aV>=0) && (aV <= gVMax) && (aD <= gDMax)) {
        if ((aD >= parabolicBound)){
            isAllowedParab = true;
        }
    }
    
    return isAllowedParab;
}

/*==============================================================================
 * Is on considered portion of the parabolic bd d_{\alpha}
 *============================================================================*/
bool isOnParabPiece (double aD, double aV, double aT){
    /*
     Function: isOnParabPiece
     
     Arguments: position (aD), velocity (aV), time (aT)
     
     Purpose: Determines if point is on parabolic portion of the boundary at time t
     
     Returns: True if on d_{\alpha}, false otherwise.
     */
    
    bool isOnParab = false;
    
    double timeToRedStart = gTR - aT; //time until red light

    double epsilon = 1e-7;
    
    double timeToRedPhase = (timeToRedStart > 0) ? timeToRedStart : 0;
    double timeToMaxV = (gVMax - aV) / gBeta; //time to reach max velocity
    
    double parabolicBound = (pow(aV,2) / (2 * gAlpha)); //max braking parab
    double maxAccelLine = (aV * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2)); //max accel line
    double maxAccelParabExtension = (aV * timeToMaxV) + (gBeta * 0.5 * pow(timeToMaxV, 2)) + (gVMax * (timeToRedPhase - timeToMaxV));
    double vCritical = gVMax - gBeta * timeToRedPhase; //velocity above which you'll reach vMax before TR if max accel
    
    //Quad form for velocity at intersection between max braking and max accel parabs
    double parabExtensionA = (gBeta + gAlpha) / (2 * gAlpha * gBeta); //a
    double parabExtensionB = -gVMax / gBeta; //b
    double parabExtensionC = pow(gVMax, 2) / (2 * gBeta) -timeToRedPhase * gVMax; //c
    double parabExtensionVelThresh = (-parabExtensionB + sqrt(pow(parabExtensionB, 2) - 4 *parabExtensionA * parabExtensionC)) / (2 * parabExtensionA); //quadratic formula

    double velocityThreshRoot = gAlpha * timeToRedPhase * (1 + sqrt(1 + (gBeta / gAlpha))); //velocity at intersection between max braking parab and max accel line
    
    double velocityThresh = (vCritical < velocityThreshRoot) ? parabExtensionVelThresh : velocityThreshRoot; //set velocity threshold
    velocityThresh = (timeToRedPhase > 0) ? velocityThresh : -epsilon;
    velocityThresh = (gFirstPass) ? -epsilon : velocityThresh; //first pass modification since no max accel line there
    
    if((fabs(aD - parabolicBound) <= epsilon) && (aV >= velocityThresh)){
        isOnParab = true;
    }
    
    return isOnParab;
}

/*==============================================================================
 * Is on considered portion of the max accel bd d_{\beta}
 *============================================================================*/
bool isOnMaxAccelBD (double aD, double aV, double aT){
    /*
     Function: isOnMaxAccelBD
     
     Arguments: position (aD), velocity (aV), time (aT)
     
     Purpose: Determines if (d,v) point falls on the max accel pw boundary portion of d_{\beta}
     
     Relies on isSubjectToMaxAccelExtension function to determine which portion of the
     piecwise boundary the point is subject to. Can only return true during the yellow phase.
     
     Returns: True if on the boundary portion of d_{\beta}, false otherwise
     */
    
    bool onMaxAccelBD = false;
    double tol = (1e-5) / gDNum;
    double timeToRedPhase = gTR - aT; //time until red phase starts
    double timeToMaxV = (gVMax - aV) / gBeta; //time until max velocity reached
    
    bool subjectToExtension = isSubjectToMaxAccelExtension(aV, aD, aT); //determine if subject to the parabolic portion of d_{\beta} or not
    
    //subject to parabolic portion of d_{\beta}
    if ((subjectToExtension) && (gHasBoundaryOn) && (aT < gTR)){
        double maxAccelParabExtension = (aV * timeToMaxV) + (gBeta * 0.5 * pow(timeToMaxV, 2)) + (gVMax * (timeToRedPhase - timeToMaxV));
        double bdDiff = fabs(maxAccelParabExtension - aD);
        
        onMaxAccelBD = (bdDiff <= tol) ? true : false;
        
    }
    
    //subject to linear portion of d_{\beta}
    else if ((!subjectToExtension) && (gHasBoundaryOn) && (aT < gTR)){
        double maxAccelLine = (aV * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2));
        double bdDiff = fabs(maxAccelLine - aD);

        onMaxAccelBD = (bdDiff <= tol) ? true : false;
        
    }
    
    return onMaxAccelBD;
}

/*========================================================================================
 * Determining if starting d,v is subject to the max accel line or the extension
 *========================================================================================*/
bool isSubjectToMaxAccelExtension(double aV, double aD, double aT){
    /*
     Function: isSubjectToMaxAccelExtension
     
     Arguments: velocity (aV), position (aD), time (aT)
     
     Purpose: Determines whether a point in the region to the left of the max accel boundary is subject to the linear or parabolic portion of d_{\beta}
     
     Returns: True if subject to parabolic portion, false otherwise
     */
    
    bool subjectToMaxAccelExtension = false;
    
    double timeToRedStart = gTR - aT; //time until red light
    double epsilon = 1e-13;
    double timeToRedPhase = (timeToRedStart > 0) ? timeToRedStart : 0;
    double timeToMaxV = (gVMax - aV) / gBeta; //time to reach max velocity
    
    double parabolicBound = (pow(aV,2) / (2 * gAlpha)); //max braking parab
    double maxAccelLine = (aV * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2)); //max accel line
    double maxAccelParabExtension = ((-pow(aV, 2)) / (2 * gBeta)) + ((aV * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeToRedPhase * gVMax); //parabolic extension when timeToMaxV < timeToRedPhase
    
    double vCritical = gVMax - gBeta * timeToRedPhase; //velocity above which you'll reach vMax via max accel before TR
    
    //Quadratic formula to solve for velocity at intersection between max accel parab and max braking parab
    double parabExtensionA = (gBeta + gAlpha) / (2 * gAlpha * gBeta); //a
    double parabExtensionB = -gVMax / gBeta; //b
    double parabExtensionC = pow(gVMax, 2) / (2 * gBeta) -timeToRedPhase * gVMax; //c
    double parabExtensionVelThresh = (-parabExtensionB + sqrt(pow(parabExtensionB, 2) - 4 *parabExtensionA * parabExtensionC)) / (2 * parabExtensionA); //quad form
    
    double velocityThreshRoot = gAlpha * timeToRedPhase * (1 + sqrt(1 + (gBeta / gAlpha))); //max accel line intersection with max braking parab
    
    double velocityThresh = (vCritical < velocityThreshRoot) ? parabExtensionVelThresh : velocityThreshRoot; //actual threshold v
    
    velocityThresh = (timeToRedPhase > 0) ? velocityThresh : -epsilon; //set velocity thresh again based on if we're in red phase or not
    
    //*******************Determining if point subject to max accel extension boundary****************************
    if (velocityThresh == parabExtensionVelThresh){
        if ((aD <= maxAccelParabExtension) && (aD >= 0) && (aV >= velocityThresh)){
            if(timeToRedPhase < timeToMaxV){
                cout << "v: " << aV << "\n";
                cout << "d: " << aD << "\n";
                cout << "time to max v: " << timeToMaxV << "\n";
                cout << "time to red: " << timeToRedPhase << "\n";
                cout << "vCritical: " << vCritical << "\n";
            }
            assert(timeToRedPhase > timeToMaxV);
            subjectToMaxAccelExtension = true;
        }
    }
    
    else if (velocityThresh == velocityThreshRoot){
        if ((aD <= maxAccelParabExtension) && (aD >= 0) && (aV > vCritical)){
            if(timeToRedPhase < timeToMaxV){
                cout << "v: " << aV << "\n";
                cout << "d: " << aD << "\n";
                cout << "time to max v: " << timeToMaxV << "\n";
                cout << "time to red: " << timeToRedPhase << "\n";
                cout << "vCritical: " << vCritical << "\n";
            }
            assert(timeToRedPhase > timeToMaxV);
            subjectToMaxAccelExtension = true;
        }
    }
    
    return subjectToMaxAccelExtension;
}

/*==============================================================================
 * Determine max accel line position for a given v
 *============================================================================*/
double maxAccelLine (double aV, double aT){
    /*
     Function: maxAccelLine
     
     Arguments: velocity (aV), time (aT)
     
     Purpose: Computes position of d_{\beta} for a given velocity
     
     Returns: Double value of location of d_{\beta}
     */
    
    double maxAccelPos = 0;
    double tol = (1e-5) / gDNum;
    double timeToRedStart = gTR - aT;
    
    //Terminal penalty time to red correction
    if(gInTerminalPenalty){
        if(aT < gTerminalT){
            double trueTimeToRed = gYellowDuration + gDT;
            double diff = fabs(trueTimeToRed - timeToRedStart);
            if(diff <= 1e-12){
                timeToRedStart = trueTimeToRed;
            }
        }
        assert(timeToRedStart <= gYellowDuration + gDT);
    }
    double timeToRedPhase = (timeToRedStart > 0) ? timeToRedStart : 0;
    double timeToMaxV = (gVMax - aV) / gBeta; //time to reach max velocity
    
    double parabolicBound = (pow(aV,2) / (2 * gAlpha)); //max braking parab
    double maxAccelLine = (aV * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2)); //max accel line
    double maxAccelParabExtension = (aV * timeToMaxV) + (gBeta * 0.5 * pow(timeToMaxV, 2)) + (gVMax * (timeToRedPhase - timeToMaxV));
    double vCritical = gVMax - gBeta * timeToRedPhase; //velocity above which you'll reach vMax before TR if max accel
    
    //Quad form for velocity at intersection between max braking and max accel parabs
    double parabExtensionA = (gBeta + gAlpha) / (2 * gAlpha * gBeta); //a
    double parabExtensionB = -gVMax / gBeta; //b
    double parabExtensionC = pow(gVMax, 2) / (2 * gBeta) -timeToRedPhase * gVMax; //c
    double parabExtensionVelThresh = (-parabExtensionB + sqrt(pow(parabExtensionB, 2) - 4 *parabExtensionA * parabExtensionC)) / (2 * parabExtensionA); //quadratic formula

    double velocityThreshRoot = gAlpha * timeToRedPhase * (1 + sqrt(1 + (gBeta / gAlpha))); //velocity at intersection between max braking parab and max accel line
    double velocityThresh = (vCritical < velocityThreshRoot) ? parabExtensionVelThresh : velocityThreshRoot; //set velocity threshold
    velocityThresh = (timeToRedPhase > 0) ? velocityThresh : -tol;
    
    if (((aV <= vCritical) && (gHasBoundaryOn) && (aT < gTR)) || ((gInTerminalPenalty) && (aV <= vCritical))){
        maxAccelPos = (aV * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2)); //max accel line
    }
    
    else if (((aV > vCritical) && (gHasBoundaryOn) && (aT < gTR)) || ((gInTerminalPenalty) && (aV > vCritical))){
        maxAccelPos = (-(pow(aV, 2)) / (2 * gBeta)) + ((aV * gVMax) / gBeta) - (0.5 * (pow(gVMax, 2) / gBeta)) + timeToRedPhase * gVMax;
    }
    
    return maxAccelPos;
}

/*==================================================================================
 * Determining velocity at intersection between max accel line and parabolic bd
 *==================================================================================*/
double threshVelocityMaxAccel(double aD, double aV, double aT){
    /*
     Function: threshVelocityMaxAccel
     
     Arguments: position (aD), velocity (aV), time (aT)
     
     Purpose: takes in d,v,t and determines the velocity at the intersection between the max accel line and parabola of last resort
     
     Returns: Double value of the velocity at the intersection between d_{\alpha} and d_{\beta}
     */
    
    double timeToRedStart = gTR - aT; //time to red phase
    double epsilon = 1e-13;
    double timeToRedPhase = (timeToRedStart > 0) ? timeToRedStart : 0;
    double timeToMaxV = (gVMax - aV) / gBeta; //time to reach max velocity
    
    double parabolicBound = (pow(aV,2) / (2 * gAlpha)); //max braking parab
    double maxAccelLine = (aV * timeToRedPhase) + (0.5 * gBeta * pow(timeToRedPhase, 2)); //max accel line
    double maxAccelParabExtension = ((-pow(aV, 2)) / (2 * gBeta)) + ((aV * gVMax) / gBeta) - (0.5 * pow(gVMax, 2) / gBeta) + (timeToRedPhase * gVMax); //parabolic extension when timeToMaxV < timeToRedPhase
    double vCritical = gVMax - gBeta * timeToRedPhase; //velocity above which you'll reach vMax before TR if max accel
    
    //Quad form for velocity at intersection between max braking and max accel parabs
    double parabExtensionA = (gBeta + gAlpha) / (2 * gAlpha * gBeta); //a
    double parabExtensionB = -gVMax / gBeta; //b
    double parabExtensionC = pow(gVMax, 2) / (2 * gBeta) -timeToRedPhase * gVMax; //c
    double parabExtensionVelThresh = (-parabExtensionB + sqrt(pow(parabExtensionB, 2) - 4 *parabExtensionA * parabExtensionC)) / (2 * parabExtensionA); //quadratic formula

    double velocityThreshRoot = gAlpha * timeToRedPhase * (1 + sqrt(1 + (gBeta / gAlpha))); //velocity at intersection between max braking parab and max accel line
    double velocityThresh = (vCritical < velocityThreshRoot) ? parabExtensionVelThresh : velocityThreshRoot; //set velocity threshold
    velocityThresh = (timeToRedPhase > 0) ? velocityThresh : -epsilon;
    
    return velocityThresh;
}

