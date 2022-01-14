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
 * File: ProbabilityFunctions.cpp
 *
 * Author: Mallory Gaspard
 *
 * Description: This file contains the functions responsible for computing
 * the conditional probability of light change at a provided time, used to
 * solve the HJB equation(s) associated with the uncertain turning yellow time
 * problems described in the Optimal Driving manuscript.
 *============================================================================*/

#include <iostream>
#include <fstream>
#include <cmath>
#include "ProbabilityFunctions.hpp"
#include "variables.h"

/*==============================================================================
 * Conditional probability determination
 *============================================================================*/
double conditionalProb(double aT){
    /*
     * Function: conditionalProb
     *
     * Arguments: time (aT)
     *
     * Purpose: computes the conditional probability of the light changing
     * at t provided the current time in the simulation
     * and an array of possible light lengths S with their associated probabilities.
     *
     * Returns a double for the probability that the light will change at t
     * provided that it has not changed before this time.
     *
     * Note: gModel2 and gPedestrianModel options are not covered in the
     * Optimal Driving manuscript.
     */
    
    double condProb = 0; //default value of light switching at time aT
    
    if (gModel1){
        condProb = discreteProb(aT);
    }
    
    else if (gModel2){
        condProb = uniformProb(aT);
    }
    
    else if (gPedestrianModel){
        condProb = gArrivalRate * gTau; //NOTE: based on current cond prob setup, need to subtract from 1
    }
    
    return condProb;
}

/*==============================================================================
 * Prob of light change for discrete light lengths
 *============================================================================*/
double discreteProb(double aT){
    /*
     Function: discreteProb
     
     Arguments: time (aT)
     
     Purpose: Computes conditional probability of light change at aT provided
     it has not changed before that time when there are discretely many
     possible turning yellow times.
     */
    
    double prob = 0;
    
    int startingIndex = -1; //index from where to start considering possible light lengths

    double probTimeT = 0; // probability of light switching at time aT
    
    double nextTimeStep = aT + gTau; // time at the next step
    assert(gTau == gDT);
    
    double tol = 1e-13;
    
    //Loop through possible light lengths to see if the current time t equals any of them
    for(int i=0; i < gNumStops; i++){
        
        double possibleLightLength = gLightTimes[i];
        double diffBetweenNextTimeAndPossibleLightLength = fabs(nextTimeStep - possibleLightLength);
        if((fabs(aT - possibleLightLength) <= tol) && (aT < possibleLightLength)){
            aT = possibleLightLength;
        }
        
        //check condition for possible light change at aT or possible light change WITHIN the next timestep
        if ((aT == possibleLightLength) || ((aT <= possibleLightLength) && (nextTimeStep > possibleLightLength) && (diffBetweenNextTimeAndPossibleLightLength >= gDT))){
            //first possible light length hit is where indexing begins
            startingIndex = i;
            probTimeT = gLightProbs[startingIndex];
            break;
        }
    }
    
    if (startingIndex >= 0){
        double rho = 0; //sum of the individual probabilities of light switching starting at time S_i

        for (int i = startingIndex; i < gNumStops; i++){
            rho += gLightProbs[i];
        }
    
        prob = probTimeT / rho;
    }
    
    else if (startingIndex == -1){
        prob = 0;
    }
    
    return prob;
}

/*==============================================================================
 * Prob of light change for uniform light length problem
 *============================================================================*/
double uniformProb(double aT){
    /*
     Function: uniformProb
     
     Arguments: time (aT)
     
     Purpose: Computes probability of light change when distribution
     on light change time is uniform.
     
     Note: This example not covered in Optimal Driving manuscript.
     */
    
    double prob = 0;
    
    double tol = gH / 1000;
    
    double currentTime = aT - gDT;
    
    prob = gDT / (gTYMax - currentTime);
    
    if (abs(prob - 1) <= tol){
        prob = 1;
    }
    
    return prob;
    
}

