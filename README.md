# Optimal Driving #

This repository contains the source code and plotting notebooks used to generate all examples presented in the Optimal Driving Under Traffic Signal Uncertainty by Mallory E. Gaspard and Alexander Vladimirsky. 


We study driver’s optimal trajectory planning under uncertainty in the duration of a traffic light’s green phase. We interpret this as an optimal control problem with an objective of minimizing the expected cost based on the fuel use, discomfort from rapid velocity changes, and time to destination. Treating this in the framework of dynamic programming, we show that the probability distribution on green phase durations gives rise to a sequence of Hamilton-Jacobi-Bellman PDEs, which are then solved numerically to obtain optimal acceleration/braking policy in feedback form. Our numerical examples illustrate the approach and highlight the role of conflicting goals and uncertainty in shaping drivers’ behavior.

# License #
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

# Manuscript #
The Optimal Driving manuscript can be found [here](https://www.sciencedirect.com/science/article/pii/S2405896322011661).

An earlier preprint version can be found [here](https://arxiv.org/abs/2201.04521).

# Contributions & Acknowledgements # 
  * The problem statement and the solution strategy were developed by Alexander Vladimirsky and Mallory E. Gaspard.
  * The implementation was carried out by Mallory E. Gaspard.
  * The manuscript was written by Mallory E. Gaspard and Alexander Vladimirsky.
  * The authors also acknowledge A. Nellis, J. van Hook, and N. Do for their preliminary work on the Phase R problem during the 2018 REU Program at Cornell University.

# Instructions #
  
## Requirements: ## 
The Optimal Driving C++ code requires the use of the following external library that must be installed in the `/usr/local/include` directory:
  * Boost - C++ library used for multidimensional array creation and storage.

The Optimal Driving results are plotted in jupyter notebooks. The plotting code requires the use of the following Python packages:
  * Numpy - Python package used for data processing.
  * Matplotlib - Python package used to produce heatmaps, contour, and trajectory plots.

## Problem Options: ##
The Optimal Driving C++ code solves the following problems which are numbered according to the examples in the Optimal Driving Under Uncertainty manuscript:

1. Indefinite green light (`stationary_green`)
2. Red-Green (`rg`)
3. Yellow-Red-Green (`yrg`)
4. Uncertain initial green duration, two possible turning green times  (`uncertain_green_2L`)
5. Uncertain initial green duration, three possible turning green times  (`uncertain_green_3L`)

## Running the Code: ##
The following instructions explain how to run the Optimal Driving Solver using the Makefile. 

To change compilers, edit `PLATFORM=` by putting your desired compiler after the `=` sign. The default compiler is set to `gcc`. 

To compile the code, type `make optDriving` in the terminal at this folder. To produce examples shown in the manuscript, type `make run EXAMPLE=(example number)`, replacing the parenthetical with the integer corresponding to the example of interest.

To delete the executable, type `make cleanMain`.

To delete any output files, type `make cleanOutput`. **Note: This is generally not recommended prior to plotting.**

# Visualizing Output # 
Results for each of the examples can be visualized by running the following jupyter notebooks: 
  * FeedbackControls.ipynb
      * Produces heatmaps of the feedback controls in the d,v domain at a specified timeslice. 
  * ContourPlots.ipynb
      * Produces contour plots of the value function in the d,v domain at a specified timeslice.  
  * DeterministicTrajectories.ipynb
      * Produces stacked plots of the optimal position, velocity, and acceleration from t = 0 to the vehicle's time of arrival at the target for Examples 1 through 3. Also produces contour plots of the feedback controls in the d,v domain at t = 0 with a marker particle indicating the vehicle's (d,v) starting point.  
  * UncertainGreenTrajectories.ipynb
      * Produces stacked plots of the optimal position, velocity, and acceleration from t = 0 to the vehicle's time of arrival at the target for Examples 4 and 5. Also produces contour plots of the feedback controls in the d,v domain at t = 0 with a marker particle indicating the vehicle's (d,v) starting point. 
  * ParetoPlotting.ipynb
      * Produces plot of the constrained pareto front between fuel cost and discomfort cost for three different maximum time to target tolerances during the indefinite green light phase. 
