# Monte Carlo program
Simple Monte Carlo program to simulate particles interacting with a Lennard-Jones potential

Compile with the fortran2008 standard for support for the execute_command_line function call
This call allows the program to create an output directory

Code was written in 2016 during my postgraduate masters course and may not comply with the current Fortran standards


## Overview of the code 

1. Get the initial coordinates of the particles and other parameters from the input files
2. Initialise any other variables needed throughout the main program e.g loop counters, summation variables and variables for physical properties.
3. The main loop of the calculation which updates the positions of the particles. Each loop a random particle is selected and moved to a new trial position. The change in energy that results from the new trial position is calculated. If this change in energy passes the Metropolis test then the trial position is accepted as the new position. If the Metropolis test fails then the trial position is rejected and the particle stays in its current position.
4. Output the data of interest e.g. the updated coordinates which could then be visualized in another program.
