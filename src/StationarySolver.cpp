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
 * File: StationarySolver.cpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains the Gauss-Jacobi and Gauss-Seidel solvers
 * used to solve the stationary HJB equation arising in the indefinite green
 * light problem. All of the examples in the manuscript utilize the Gauss-Seidel
 * solver with the Modified Gauss-Seidel option set to `true`. The Modified
 * Gauss-Seidel approach allows us to decouple the system of equations to
 * converge faster.
 *
 * Both methods begin with the initial data: u(d,v) = (d - gDTarget)^2
 * which is set in the greenLightOnly function in ProblemFunctions.cpp
 *============================================================================*/

//----------------Libraries--------------------------------------------------
#include "StationarySolver.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include "boost/multi_array.hpp"

//-------------------Problem specific header files----------------------------
#include "variables.h"
#include "InitFunctions.hpp"
#include "CoordinateFunctions.hpp"
#include "ProbabilityFunctions.hpp"
#include "HelperFunctions.hpp"

/*==============================================================================
 * Stationary Green Light Solver - Gauss Jacobi
 *============================================================================*/
void stationaryGJSolver(multiarray *aU, multiarray *aOC){
    /*
     Function: stationaryGJSolver
     
     Arguments: multiarrays to hold value function and optimal controls
     
     Purpose: Solves the stationary green light problem using Gauss Jacobi iteration.
     The initial guess (d - gDTarget)^2 is supplied in the greenLightOnly function in the
     ProblemFunctions.cpp file.
     
     Stencil:
     
     uTopLeft----------------uTopRight
        |                        |
        |                        |
        |                        |
        |                        |
        |                        |
        |                        |
        |                        |
        |                        |
     uBottomLeft-------------uBottomRight
     
     */
    
    gInStationarySolve = true;
    
    cout << "In stationary green GJ solver. " << "\n";
    
    double convTol = 1e-5; //convergence tolerance
    double vfDiff = gInfty; //starting difference
    short int k = 0; //iteration counter
    
    //while loop until convergence
    while (vfDiff > convTol){
        for(short int i = 1; i < gDNum + 1; i++) {     // distance loop
            double dCurrent = i * gDPos + gDTarget; //physical pos
            for(short int j = gVNum; j >= 0; j--) {   // vel loop - sweep down
                double vCurrent = j * gDVel; //physical vel
                
                //Loop over controls to find min
                double optimalControl;
                double bestValueFunction = gInfty;
                
                for (int controlIndex = 0; controlIndex < gNumControls; controlIndex++){
                    
                    double currentValue; //value function value placeholder
                    double currentControl = gControlValues[controlIndex];
                    double tau = gDT;
                    
                    //Compute dNext, vNext
                    double dNext = dCurrent - (vCurrent * tau + 0.5 * currentControl * pow(tau, 2));
                    double vNext = vCurrent + currentControl * tau;
                    
                    double dRescaled = (dNext - gDTarget) / gDPos;
                    double vRescaled = vNext / gDVel;
                    
                    double dC = ceil(dRescaled);
                    double dF = dC - 1;
                    double vF = floor(vRescaled);
                    double vC = vF + 1;
                    
                    int dCIdx = int(dC);
                    int dFIdx = dCIdx - 1;
                    int vFIdx = int(vF);
                    int vCIdx = vFIdx + 1;
                    
                    //adjustment if j == gvMax -- downward bias
                    vC = (j == gVNum) ? j : vC;
                    vF = (j == gVNum) ? j - 1 : vF;
                    vCIdx = (j == gVNum) ? j : vCIdx;
                    vFIdx = (j == gVNum) ? j - 1 : vFIdx;
                    
                    //Assign interp weights - bilinear interp
                    double alphaBottomRight = ((dRescaled - dF) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
                    double alphaBottomLeft = ((dC - dRescaled) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
                    double alphaTopLeft = ((dC - dRescaled) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));
                    double alphaTopRight = ((dRescaled - dF) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));
                    
                    double alphaSum = alphaTopRight + alphaBottomRight + alphaTopLeft + alphaBottomLeft;
                    assert(fabs(alphaSum - 1) <= 1e-12);
                    
                    if((vNext > gVMax) || (vNext < 0) || (dNext > gDMax) || (dNext < gDTarget)){
                        currentValue = gInfty;
                    }
                    
                    else{
                        double uBottomLeft = (*aU)[dFIdx][vFIdx][0];
                        double uTopLeft = (*aU)[dFIdx][vCIdx][0];
                        double uBottomRight = (*aU)[dCIdx][vFIdx][0];
                        double uTopRight = (*aU)[dCIdx][vCIdx][0];
                        
                        double runningCost = runningCostIntegral(currentControl, 0, tau, gGamma, gC1, gC2, gC3);
                        
                        if(j == gVNum){
                            currentValue = runningCost + alphaTopLeft * uTopLeft + alphaBottomLeft * uBottomLeft + alphaBottomRight * uBottomRight + alphaTopRight * uTopRight;
                        }
                        
                        else{
                            currentValue = runningCost + alphaTopLeft * uTopLeft + alphaBottomLeft * uBottomLeft + alphaTopRight * uTopRight + alphaBottomRight * uBottomRight;
                        }
                    }//end else
                    
                    //Comapre values
                    if(currentValue < bestValueFunction){
                        bestValueFunction = currentValue;
                        optimalControl = currentControl;
                        assert(bestValueFunction < gInfty);
                        assert(bestValueFunction >= 0);
                    }
                    
                }//end loop over controls
                
                //Save new value to the 1-index spot
                (*aU)[i][j][1] = bestValueFunction;
                (*aOC)[i][j][1] = optimalControl;
                assert(optimalControl >= 0);
                
            }//end loop over j
        }//end loop over i
        
        double infNorm = infinityNorm(aU, -1); //compute infinity norm to to check for convergence
        
        vfDiff = infNorm;
        cout << "current norm: " << vfDiff << "k:" << k << "\n";

        //reset values for next iteration
        for(short int i = 0; i < gDNum + 1; i++) {     // distance loop
            for(short int j = 0; j < gVNum + 1; j++){ // velocity loop
                double newVFVal = (*aU)[i][j][1];
                double newOCVal = (*aOC)[i][j][1];
                
                (*aU)[i][j][0] = newVFVal;
                (*aOC)[i][j][0] = newOCVal;
            }
        }
        
        k++;
    }//end while loop over k
    cout << "Through stationary solve" << "\n";
    gInGJ = false;
}

/*==============================================================================
 * Stationary Green Light Solver - Gauss Seidel
 *============================================================================*/
void stationaryGSSolver(multiarray *aU, multiarray *aOC){
    /*
     Function: stationaryGSSolver()
     
     Purpose: Solves the stationary green light problem using gauss seidel iteration.
     Has option for Modified GS or regular GS although all of the examples shown in the
     manuscript use Modified GS. The Modified GS method eliminates the RHS dependence
     on the current node's Q_{i,j} value which significantly speeds up convergence.
     
     This function relies on the goldenSection function for determining the optimal
     control and value function via GSS after the control grid-search, and the infinityNorm
     function to compute the infinity norm over all j's for a fixed i index to test for
     convergence.
     
     Stencil:
     
     uTopLeft----------------uTopRight
        |                        |
        |                        |
        |                        |
        |                        |
        |                        |
        |                        |
        |                        |
        |                        |
     uBottomLeft-------------uBottomRight
     
     */
    
    gInStationarySolve = true;
    gInGS = true;
    
    cout << "In Gauss Seidel." << "\n";
    
    double convTol = 1e-5; //convergence tolerance
    short int dMaxIdx = (gSolvingGreenOnly) ? gDNum + 1 : (gDNum / 2) + 1;
    
    for(short int i = 1; i < dMaxIdx; i++){
        double dCurrent = (gSolvingGreenOnly) ? i * gDPos + gDTarget : i * gDPos;
        gPosIndex = i;
        short int k = 0; //iteration index
        double maxNorm = gInfty;
        
        while (maxNorm > convTol){
            gTimeIndex = k;
            for(short int j = gVNum; j >= 0; j--){
                gVelIndex = j;
                double vCurrent = j * gDVel;
                
                double bestValueFunction = gInfty;
                double optimalControl = 0.0;

                if(gSolvingGreenOnly){
                    for (int controlIndex = 0; controlIndex < gNumControls; controlIndex++){
                        
                        double currentValue; //value function value placeholder
                        double currentControl = gControlValues[controlIndex];
                        double tau = gDT;
                        
                        //Compute dNext, vNext
                        double dNext = dCurrent - (vCurrent * tau + 0.5 * currentControl * pow(tau, 2));
                        double vNext = vCurrent + currentControl * tau;
                        
                        double dRescaled = (dNext - gDTarget) / gDPos;
                        double vRescaled = vNext / gDVel;
                        
                        double dC = i;
                        double dF = dC - 1;
                        double vF = j;
                        double vC = vF + 1;
                        
                        int dCIdx = int(dC);
                        int dFIdx = dCIdx - 1;
                        int vFIdx = int(vF);
                        int vCIdx = vFIdx + 1;
                        
                        //adjustment if j == gvMax -- downward bias
                        vC = (j == gVNum) ? j : vC;
                        vF = (j == gVNum) ? j - 1 : vF;
                        vCIdx = (j == gVNum) ? j : vCIdx;
                        vFIdx = (j == gVNum) ? j - 1 : vFIdx;
                        
                        //Assign interp weights

                        double alphaBottomRight = ((dRescaled - dF) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
                        double alphaBottomLeft = ((dC - dRescaled) * (vC - vRescaled)) / ((dC - dF) * (vC - vF));
                        double alphaTopLeft = ((dC - dRescaled) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));
                        double alphaTopRight = ((dRescaled - dF) * (vRescaled - vF)) / ((dC - dF) * (vC - vF));
                        
                        alphaTopRight = (fabs(alphaTopRight) <= 1e-13) ? 0 : alphaTopRight;
                        alphaTopLeft = (fabs(alphaTopLeft) <= 1e-13) ? 0 : alphaTopLeft;
                        alphaBottomRight = (fabs(alphaBottomRight) <= 1e-13) ? 0 : alphaBottomRight;
                        alphaBottomLeft = (fabs(alphaBottomLeft) <= 1e-13) ? 0 : alphaBottomLeft;
                        
                        alphaTopRight = (fabs(alphaTopRight - 1) <= 1e-13) ? 1 : alphaTopRight;
                        alphaTopLeft = (fabs(alphaTopLeft - 1) <= 1e-13) ? 1 : alphaTopLeft;
                        alphaBottomRight = (fabs(alphaBottomRight - 1) <= 1e-13) ? 1 : alphaBottomRight;
                        alphaBottomLeft = (fabs(alphaBottomLeft - 1) <= 1e-13) ? 1 : alphaBottomLeft;
                        
                        double alphaSum = alphaTopRight + alphaBottomRight + alphaTopLeft + alphaBottomLeft;
                        assert(fabs(alphaSum - 1) <= 1e-12);
                        
                        if((vNext > gVMax) || (vNext < 0) || (dNext > gDMax) || (dNext < gDTarget)){
                            currentValue = gInfty;
                        }
                        
                        if(j == gVNum){
                            //NOTE: Hard coded horizontal linear interp for a^* = 0 along this line!
                            vC = j;
                            vF = vC;
                            vCIdx = j;
                            vFIdx = vCIdx;
                            
                            alphaTopLeft = (dRescaled - dF) / (dC - dF);
                            alphaTopRight = 1 - alphaTopLeft;
                            
                            double uTopLeft = (*aU)[i-1][j][0]; //converged value on i-1
                            double uTopRight = (*aU)[i][j][0]; //current guess
                            
                            double runningCost = runningCostIntegral(0, 0, tau, gGamma, gC1, gC2, gC3); //Note: current time arg is irrelevant

                            currentValue = (1 / (1 - alphaTopLeft)) * (runningCost + (1 - alphaTopLeft) * uTopLeft);
                            
                            bestValueFunction = currentValue;
                            optimalControl = 0; //0 must be OC here, set manually

                            if(k >= 1){
                                assert(uTopRight == currentValue);
                            }
                            
                            break;
                        }
                        
                        else{
                            //Left values have converged - take from 0 since these columns have been written with the converged solution
                            double uBottomLeft = (*aU)[dFIdx][vFIdx][0];
                            double uTopLeft = (*aU)[dFIdx][vCIdx][0];
                            
                            double uBottomRight = (*aU)[dCIdx][vFIdx][1]; //current guess from before
                            double uTopRight = (*aU)[dCIdx][vCIdx][1]; //previous iteration update at j+1 saved in 1 spot
                            
                            double runningCost = runningCostIntegral(currentControl, 0, tau, gGamma, gC1, gC2, gC3); //Note: current time arg is irrelevant
                            
                            if (gModifiedGS){
                                double gamma = (1 / (1 - alphaBottomRight));
                                currentValue = gamma * (runningCost + alphaTopLeft * uTopLeft + alphaTopRight * uTopRight + alphaBottomLeft * uBottomLeft);
                                
                            }
                            
                            else{
                                currentValue = runningCost + alphaTopLeft * uTopLeft + alphaBottomLeft * uBottomLeft + alphaTopRight * uTopRight + alphaBottomRight * uBottomRight;
                            }
                            
                        }//end else
                        
                        //Comapre values
                        if(currentValue < bestValueFunction){
                            bestValueFunction = currentValue;
                            optimalControl = currentControl;
                            assert(bestValueFunction < gInfty);
                            assert(bestValueFunction >= 0);
                        }
                        
                    }//end loop over controls
                    
                    //Call GSS
                    if((j!= gVNum) && (i >= 1)){
                    
                        double searchIntervalWidth = gControlSpacing;
                        double aLeft = (optimalControl - searchIntervalWidth < -gAlpha) ? -gAlpha : optimalControl - searchIntervalWidth;
                        if((gVelIndex == 0) || (optimalControl == 0)){
                            aLeft = 0;
                        }

                        double aRight = (optimalControl + searchIntervalWidth > gBeta) ? gBeta : optimalControl + searchIntervalWidth;
                        
                        assert(aLeft >= -gAlpha);
                        assert(aRight <= gBeta);
                        
                        struct sGSS gssValues = goldenSection(aU, aOC, aU, aOC, aU, aOC, 0, aLeft, aRight, 0, dCurrent, vCurrent, gDT);
                        
                        bestValueFunction = gssValues.valueFn;
                        optimalControl = gssValues.oc;
                        
                    }
                    
                    (*aU)[i][j][1] = bestValueFunction;
                    (*aOC)[i][j][1] = optimalControl;
                }//end solving green only
            }//end loop over j
            
            //Compute infinity norm over j
            double infNorm = infinityNorm(aU, i); //get infinity norm of i^th row
            
            maxNorm = infNorm;

            //Replace i^th row of aU and aOC with the converged values over j's
            for(short int j = 0; j < gVNum + 1; j++){
                double newVal = (*aU)[i][j][1];
                assert(newVal < gInfty);
                assert(newVal >= 0);
                (*aU)[i][j][0] = (*aU)[i][j][1];
                if(gSolvingGreenOnly){
                    (*aOC)[i][j][0] = (*aOC)[i][j][1];
                }
            }
            
            k++;
        }//end while loop
        
    }//end loop over i
    
    cout << "Out of Gauss Seidel." << "\n";
    gInGS = false;
}


