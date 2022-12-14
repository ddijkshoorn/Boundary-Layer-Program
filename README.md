# Boundary-Layer-Program
## Computer Program for the Simulation of Two-dimensional Steady State Boundary Layers in Non-ideal Gas Flows

Dominic Dyon Dijkshoorn, 14-12-2022

This MATLAB computer program solves the two-dimensional steady state boundary layer equations with general fluid properties for compressible flows in the ideal gas or non-ideal gas (departing from ideal gas) regime, adiabatic or including heat transfer for laminar and/or turbulent (algebraic CS-model) flows. The computer program was verified and validated for air, and needs validation for flows departing from ideal gas.

The core of the computer program (method of solution) was taken from one of the FORTRAN programs contained on the DVD enclosed with the book Convective Heat Transfer (2002) by Tuncer Cebeci [3]. The FORTRAN code has been converted to MATLAB code, and improved by adding an omitted coefficient. The simple computer program utilizing calorically perfect ideal gas was extended to include three gas models: calorically perfect ideal gas with constant fluid properties (but with viscosity as function of temperature), thermally perfect (calorically imperfect) ideal gas with fluid properties as function of temperature, and, non-ideal gas retrieving fluid properties with state-of-the-art thermophysical models through FluidProp [5]. The pre-determined point-of-transition method for transition from laminar to turbulent flow was extended with two simple transition prediction methods for air (Wazzan’s and Michel’s method), and a re-laminarisation prediction method (Nash-Webber), also for air (see thesis [1] for more details on these methods). The algebraic turbulence model implemented (CS-model) was extended with a turbulent Prandtl-number model obtained from the book Analysis of Turbulent Boundary Layers (1974) by Cebeci and Smith [4].

This MATLAB computer program was developed and used by Dijkshoorn [1] to numerically investigate the formation of a boundary layer along a wall surface inside a converging-diverging nozzle accelerating Hexamethyldisiloxane (MM) departing from ideal gas (compressibility factor Z=0.55 at maximum departure from ideal gas) for the ORCHID test set-up [2] as part of the master thesis Simulation of Two-Dimensional Steady State Boundary Layers Applied to Nonideal Gas Flows, which was completed on August 10, 2020 at Delft University of Technology, faculty of Mechanical, Maritime and Materials Engineering (3me) [1]. The thesis includes an extensive explanation of the computer program, including theory and the structure of the program itself, verification and validation for air, and the aforementioned numerical investigation of the non-ideal gas expansion of MM inside a converging-diverging nozzle. The thesis can be found on the TU-Delft repository, via the following link: https://repository.tudelft.nl/islandora/object/uuid%3A433d4c00-e063-4614-9c78-2849870d8d7f

## Structure of the program and file-folders
The computer program itself, consisting of different function files, is stored in the main-folder. Thesis figure 3-2 [1] visualises the structure of the program using a flow-chart. The flow-chart relates all function files in the main-folder. Table 1 in the Addendum [1] gives an overview of the variables and parameters and accompanying structures. Auxiliary files, such as input-files and data-files, are stored in subfolders named accordingly. The verification and validation process and an accuracy study are also included in subfolders named accordingly. All files can be re-run to obtain the data presented in the thesis. The data files resulting from the simulations are stored for completeness in a subfolder accordingly (folder ‘Stored_Sim_Data’).

## Running an example case (NACA0012 airfoil)
Open and run the MAIN-file in the main folder. Graphs will be generated to present the simulation results. The example case, a boundary layer flow over a NACA0012 airfoil, is obtained with the FORTRAN program from the aforementioned DVD [3], since it was integrated in the example FORTAN code. The example treats incompressible stagnation point flow following the top-side of a NACA0012 airfoil (Ma<<0.3), including laminar flow, transition to turbulent flow, turbulent flow, and lastly, flow separation. To run a different case: replace the input-file in the MAIN-file by any other input-file from the 'INPUT'-folder or 'Validation_Study/INPUT'-folder. The folder ‘Simulation’ contains the caller-files for all non-ideal gas simulations from chapter 5 in the master thesis.

NB for the non-ideal gas simulations the licensed program FluidProp [5] is required to run the simulations. After having acquired FluidProp, copy the following files to the MAIN-folder: ‘InitFluidProp.m’ and ‘Cleanup_FluidProp.m’. 

## References
[1] Dijkshoorn, D.D. (2020). Simulation of Two-Dimensional Steady State Boundary Layers Applied to Nonideal Gas Flows. Master’s Thesis, Delft University of Technology, faculty of Mechanical, Maritime and Materials Engineering (3me), The Netherlands. Link:
https://repository.tudelft.nl/islandora/object/uuid%3A433d4c00-e063-4614-9c78-2849870d8d7f 

[2] Head, A. J. (2021). Novel experiments for the investigation of non-ideal compressible fluid dynamics: the ORCHID and first results of optical measurements. PhD Thesis, Faculty of Aerospace Engineering, Delft University of Technology, The Netherlands.

[3] Cebeci, T. (2002). Convective heat transfer. 2nd ed. Springer Berlin Heidelberg & Horizons Publishing, Long Beach California.

[4] Cebeci, T. and Smith, A.M.O. (1974). Analysis of turbulent boundary layers. Academic Press, New York.

[5] Colonna, P. and Van der Stelt T.P. (2004). FluidProp: a program for the estimation of thermo physical properties of fluids. Energy Technology Section, Delft University of Technology, The Netherlands.
