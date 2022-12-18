# Boundary-Layer-Program

## Computer Program for the Simulation of Two-dimensional Steady State Boundary Layers in Non-ideal Gas Flows

*Dominic Dyon Dijkshoorn, 14-12-2022*

This MATLAB computer program solves the two-dimensional steady state boundary layer equations with general fluid properties for compressible flows in the ideal gas or non-ideal gas (departing from ideal gas) regime, adiabatic or including heat transfer for laminar and/or turbulent (algebraic CS-model) flows. The computer program was verified and validated for air, and needs validation for flows departing from ideal gas.

The core of the computer program (method of solution) was taken from one of the FORTRAN programs contained on the DVD enclosed with the book *Convective Heat Transfer* (2002) by Tuncer Cebeci [3]. The FORTRAN code has been converted to MATLAB code, and improved by adding an omitted coefficient. The simple computer program utilizing calorically perfect ideal gas was extended to include three gas models: calorically perfect ideal gas with constant fluid properties (but with viscosity as function of temperature), thermally perfect (calorically imperfect) ideal gas with fluid properties as function of temperature, and, non-ideal gas retrieving fluid properties with state-of-the-art thermophysical models through FluidProp [5]. The pre-determined point-of-transition method for transition from laminar to turbulent flow was extended with two simple transition prediction methods for air (Wazzan’s and Michel’s method), and a re-laminarisation prediction method (Nash-Webber), also for air (see thesis [1] for more details on these methods). The algebraic turbulence model implemented (CS-model) was extended with a turbulent Prandtl-number model obtained from the book *Analysis of Turbulent Boundary Layers* (1974) by Cebeci and Smith [4].

This MATLAB computer program was developed and used by Dijkshoorn [1] to numerically investigate the formation of a boundary layer along a wall surface inside a converging-diverging nozzle accelerating Hexamethyldisiloxane (MM) departing from ideal gas (compressibility factor *Z*=0.55 at maximum departure from ideal gas) for the ORCHID test set-up [2] as part of the master thesis *Simulation of Two-Dimensional Steady State Boundary Layers Applied to Nonideal Gas Flows*, which was completed on August 10, 2020 at Delft University of Technology, faculty of Mechanical, Maritime and Materials Engineering (3me) [1]. The thesis includes an extensive explanation of the computer program, including theory and the structure of the program itself, verification and validation for air, and the aforementioned numerical investigation of the non-ideal gas expansion of MM inside a converging-diverging nozzle. The thesis can be found on the TU-Delft repository, via the following link: https://repository.tudelft.nl/islandora/object/uuid%3A433d4c00-e063-4614-9c78-2849870d8d7f

## Structure of the program and file-folders

The computer program itself, consisting of different function files, is stored in the main-folder. Thesis figure 3-2 [1] visualises the structure of the program using a flow-chart. The flow-chart relates all function files in the main-folder. Table 1 below (end of page below references, also included in the Addendum [1]) gives an overview of the variables and parameters and accompanying structures. Auxiliary files, such as input-files and data-files, are stored in subfolders named accordingly. The verification and validation process and an accuracy study are also included in subfolders named accordingly. All files can be re-run to obtain the data presented in the thesis. The data files resulting from the simulations are stored for completeness in a subfolder accordingly (folder ‘Stored_Sim_Data’).

## Running an example case (NACA0012 airfoil)

Open and run the MAIN-file in the main folder. Graphs will be generated to present the simulation results. The example case, a boundary layer flow over a NACA0012 airfoil, is obtained with the FORTRAN program from the aforementioned DVD [3], since it was integrated in the example FORTAN code. The example treats incompressible stagnation point flow following the top-side of a NACA0012 airfoil (Ma<<0.3), including laminar flow, transition to turbulent flow, turbulent flow, and lastly, flow separation. The To run a different case: replace the input-file in the MAIN-file by any other input-file from the 'INPUT'-folder or 'Validation_Study/INPUT'-folder. The folder ‘Simulation’ contains the caller-files for all non-ideal gas simulations from chapter 5 in the master thesis.

NB for the non-ideal gas simulations the licensed program FluidProp [5] is required to run the simulations. After having acquired FluidProp, copy the following files to the MAIN-folder: ‘InitFluidProp.m’ and ‘Cleanup_FluidProp.m’. 

## References

[1] Dijkshoorn, D.D. (2020). *Simulation of Two-Dimensional Steady State Boundary Layers Applied to Nonideal Gas Flows.* Master’s Thesis, Delft University of Technology, faculty of Mechanical, Maritime and Materials Engineering (3me), The Netherlands. Link:
https://repository.tudelft.nl/islandora/object/uuid%3A433d4c00-e063-4614-9c78-2849870d8d7f 

[2] Head, A. J. (2021). *Novel experiments for the investigation of non-ideal compressible fluid dynamics: the ORCHID and first results of optical measurements.* PhD Thesis, Faculty of Aerospace Engineering, Delft University of Technology, The Netherlands.

[3] Cebeci, T. (2002). *Convective heat transfer.* 2nd ed. Springer Berlin Heidelberg & Horizons Publishing, Long Beach California.

[4] Cebeci, T. and Smith, A.M.O. (1974). *Analysis of turbulent boundary layers.* Academic Press, New York.

[5] Colonna, P. and Van der Stelt T.P. (2004). *FluidProp: a program for the estimation of thermo physical properties of fluids.* Energy Technology Section, Delft University of Technology, The Netherlands.

## Table

*Table 1: List of structures and variables used in the MATLAB program. NB The calculation method of the most important boundary layer characteristics has been kept the same for the purpose of verification with the FORTRAN program. All variables include comments in code.*
| Structure   | Variable           | Description                                                                                                |
| :---        | :---               | :---                                                                                                       |
| BLC         |                    | BL Characteristics: contains all calculated Boundary Layer Characteristics (properties, vector or matrix)  |
|             | edge               | BL Edge location\height (vector)                                                                           |
|             | delta              | BL velocity thickness                                                                                      |
|             | delta_ast          | BL displacement thickness                                                                                  |
|             | theta              | BL momentum thickness                                                                                      |
|             | delta3             | BL kinetic energy thickness                                                                                |
|             | delta4             | BL enthalpy thickness (according to different definitions, see comments in code)                           |
|             | delta5             | BL entropy thickness                                                                                       |
|             | H                  | Boundary Layer Shape factor H                                                                              |
|             | Re_x               | Reynolds-number                                                                                            |
|             | Re_delta_ast       | Reynolds-number displacement thickness                                                                     |
|             | Re_theta           | Reynolds-number momentum thickness                                                                         |
|             | Pr_x               | Prandtl-number                                                                                             |
|             | Pe_x               | Peclet-number                                                                                              |
|             | St_x               | Stanton-number                                                                                             |
|             | Nu_x               | Nusselt-number                                                                                             |
|             | Ec_x               | Eckert-number                                                                                              |
|             | Br_x               | Brinkman-number                                                                                            | 
|             | Cd                 | Denton's loss coefficient                                                                                  |
|             |                | 
|             |                | 
|             |                | 
|             |                | 
|             |                | 
| EDG         |                    | BL edge: contains all Boundary Layer Edge properties (every varibale is a vector with values along entire BL Edge) |
|             | HtE                | Total Enthalpy at BL Edge                                                                                  |
|             | TtE                | Total Temperature BL Edge                                                                                  |
|             | PtE                | Total Pressure at BL Edge                                                                                  |
|             | TsE                | Static Temperature BL Edge                                                                                 |
|             | UE                 | Flow velocity at the Boundary Layer Edge                                                                   |
|             | CpE                | Cp at BL Edge                                                                                              |
|             | HsE                | Static Enthalpy at BL Edge                                                                                 |
|             | muE                | Viscosity at BL Edge                                                                                       |
|             | kE                 | Conductivity at BL Edge                                                                                    |
|             | gammaE             | gamma at BL Edge                                                                                           |
|             | aE                 | Speed of sound (SoS) at BL Edge                                                                            |
|             | MaE                | Mach-number at BL Edge                                                                                     |
|             | PsE                | Static Pressure at BL Edge                                                                                 |
|             | rhoE               | Static density at BL Edge                                                                                  |
|             | Re_x               | Reynolds-number along BL Edge                                                                              |
| FLP         |                    | FLuid Properties: contains all FLuid Properties inside the BL                                              |
|             |                 |                                                                                |
|             |                 |                                                                                |
|             |                 |                                                                                |
|             |                 |                                                                                |
|             |                 |                                                                                |
|             |                 |                                                                                |
|             |                 |                                                                                |
| FRS         |                    | FRee Stream Initial conditions: contains all FRee Stream properties (outside BL) at the start              |
|             | PsI                | Initial Static Pressure                                                                                    |
|             | TsI                | Initial Static Temperature                                                                                 |
|             | UI                 | Initial flow velocity                                                                                      |
|             | gammaI             | gamma for above flow conditions                                                                            |
|             | CpI                | Cp for above flow conditions                                                                               |
|             | HsI                | Initial static enthalpy                                                                                    |
|             | MaI                | Initial Mach-number                                                                                        |
|             | aI                 | Initial speed of sound                                                                                     |
|             | TtI                | Initial total Temperature                                                                                  |
|             | HtI                | Initial total Enthalpy                                                                                     |
|             | PtI                | Initial total Pressure                                                                                     |
|             | rhoI               | Initial density (static)                                                                                   |
|             | muI                | Initial viscosity                                                                                          |
|             | kI                 | Initial conductivity                                                                                       |
|             | Re_I               | Initial Reynolds-number                                                                                    |
|             | Pr_I               | Initial Prandtl-number                                                                                     |
|             | Pe_I               | Initial Peclet-number                                                                                      |
|             | Ec_I               | Initial Eckert-number (for comparison)                                                                     |
|             | Br_I               | Initial Brinkman-number                                                                                    |
| GRD         |                    | GRiD: contains all GRiD properties                                                                         |
|             | etaE               | Maximum grid height (grid height at BL Edge)                                                               |
|             | VGP                | Variable Grid Parameter: multiplication factor which determines the grid spacing                           |
|             | Deta               | Differences between the eta-grid points (vertical grid): all differences in vector form                    |
|             | NP                 | Number of grid-Points in vertical direction                                                                |
|             | eta                | The eta-grid (vector)                                                                                      |
|             | A                  | Help variable in solution method                                                                           |
| HVR         |                    | Help VaRiables: contains all property derived variables for solving the set of differential equations      |
|             | P1                 | Pressure Gradient Parameter m1                                                                             |
|             | P2                 | Pressure Gradient Parameter m2                                                                             |
|             | CEL                | Help variable representing dimensionless *X*-coordinate for solver                                           |
|             | P1P                | Pressure Gradient Parameter adjusted with CEL to solver coefficient                                        |
|             | P2P                | Pressure Gradient Parameter adjusted with CEL to solver coefficient                                        |
|             | alpha0             | Solver coefficient determining adiabatic or heat transfer (0 or 1)                                         |
|             | alpha1             | Solver coefficient determining adiabatic or heat transfer (0 or 1)                                         |
|             | WW                 | Boundary condition solver input                                                                            |
| INP         |                    | INPut: contains all standard and case specific INPut properties; combination dependent on sepcific case: 3 initial free stream conditions properties (flat plate flow or stagnation point flow) and 1 EDG input variable |
|             | x                  | Dimensionless *x*-coordinate of surface geometry                                                           |
|             | y                  | Dimensionless *y*-coordinate of surface geometry                                                           |
|             | L                  | Scale factor (not used, equal to 1) for scaling *x* and *y*                                                |
|             | PsI                | Initial static Pressure                                                                                    |
|             | TsI                | Initial static Temperature                                                                                 |
|             | UI                 | Initial velocity                                                                                           |
|             | uE                 | BL Edge velocity (vector)                                                                                  |
|             | BCW                | Boundary Condition Wall                                                                                    |
| MON         |                    | MONitor: contains all variables that count or monitor the numerical process properties                     |
|             | ITE                | Number of ITErations per station                                                                           |
|             | MGE                | Total number of Mesh Grid Extensions per station                                                           |
|             | NST                | Total Number of STations                                                                                   |
|             | NP                 | Amount of vertical grid points per station (vector); this number grows with mesh grid extensions           |
|             | SEP                | SEParation has taken place yes/no                                                                          |
|             | STR                | Station at which Separation has taken place                                                                |
|             | NTR                | StatioN at which TRansition has taken place                                                                |
|             | Re_tr              | Reynolds-number at station where transition took place                                                     |
|             | tr                 | Transition (tr) has taken place yes/no                                                                     |
|             | tr_last_check      | Station at which last transition check took place                                                          |
|             | rl                 | Re-laminarisation (rl) has taken place yes/no                                                              |
|             | rl_last_check      | Station at which last re-laminarisation check took place                                                   |
| OPT         |                    | OPTions: contains all input OPTions available defining the case; see INPUT-file for a detailed explanation |
|             | GASM               | Gas model: calorically perfect IG; Thermally perfect IG; or, Nonideal Gas                                  |
|             | COMP               | Compressible or incompressible flow                                                                        |
|             | CPRN               | Constant or variable Prandtl-number                                                                        |
|             | CCRP               | Constant or variable/general Chapman-Rubesin parameter                                                     |
|             | CPRT               | Constant or variable turbulent Pr-number                                                                   |
|             | BCIE               | Boundary Condition Input Boundary Layer Edge (BL Edge input options: 1=uE (UE/UI); 2=MaE; 3=psE (PsE/PtI)) |
|             | BCEE               | Wall Boundary Condition Energy Equation: adiabatic wall or heat transfer                                   |
|             | TRME               | Transition method: no transition; prescribed transition location NTR; Wazzan's method; Michel's method; or, fully turbulent    |
|             | RLAM               | Re-laminarization method, no relam.; prescribed location NRL; or, simple engineering estimate based on exp. data for air only! |
|             | GRAD               | Method to calculate derivative of a vector variable (Lagrange; Weighted-difference technique; or, SPLINE methods)              |
| PLT         |                    | PLoT: contains all PLoTting parameters of what to plot                                                     |
| SET         |                    | SETtings: contains all numerical settings; see INPUT-file for a detailed explanation                       |
|             | NPT                | Maximum number of grid-points (stations) in *X*-direction                                                  |
|             | NTR                | Pre-defined transition station                                                                             |
|             | NRL                | Pre-defined laminarisation station                                                                         |
|             | ITMAX              | Maximum number of iterations in main calculation loop                                                      |
|             | ITMAX0             | Maximum number of iterations for initial profile calculation in IVPL-file                                  |
|             | ITMAX2             | Maximum number of iterations in gas model 2 fluid property calculation (Newton-Raphson) in PRECAL-file     |
|             | ITMAX3             | Maximum number of iterations in gas model 3 fluid property calculation (Newton-Raphson) in PRECAL-file     |
| SOL         |                    | SOLution: contains all SOLution parameters (solver variables: f, u, v, g, p, b, c, d, e)                   |
| TCC         |                    | Turbulent CharaCteristics: contains all parameters and constants (closure coefficients) involved in the algebraic turbulence model |
|             | kappa              | Von Karman constant                                                                                        |
|             | A_plus             | Van Driest damping constant/factor                                                                         |
|             | alpha              | Clauser's/outer eddy viscosity constant                                                                    |
|             | ints               | Assumed intersection of viscous sublayer with (intermediate) log layer (Cebeci, 1974)                      |
|             | kappa_h            | Heat transfer mixing-length constant                                                                       |
|             | PrT                | Turbulent Prandtl-number (constant or vector)                                                              |
|             | Bcoeff             | Constants in fluid specific eddy conductivity damping factor (Turbulent Prandtl-number model)              |
|             | gamma_tr           | Transition region intermittency factor (for laminar to turbulent flow transition)                          |
|             | gamma_int          | Klebanoff intermittency factor                                                                             |
