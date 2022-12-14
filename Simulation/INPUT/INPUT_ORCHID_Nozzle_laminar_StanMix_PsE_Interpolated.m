%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT File                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Written by D.D.Dijkshoorn (DDD) on 31/01/2020             %
%   Last checked: 09/03/2022                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Case: de Laval nozzle designed for ORCHID test section
%   Source: Adam Head (design and SU2 simulation results (euler))
%   from surface-flow Excel-file PR0250_EULER_StanMix_15K_2nd_20200227 (10-03-2020)
%   Verification of: transsonic nozzle flow simulations
%   Goal: estimate Boundary Layer characteristics of ORCHID Nozzle design:
%   Do we need to take into account the viscous effects by adding the
%   displacement thickness?

% Simulation StanMix (SM) takes about 2 minutes

Case_name = 'ORCHID_Nozzle_StanMix_PsE_Interpolated_lam'; % tag

%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ideal Gas (IG); Nonideal Gas (NG), Real Gas (RG) (not applied here)

OPT.GASM = 3;       % Gas model: 1=calorically perfect IG; 2=thermally perfect IG; 3=Nonideal Gas (NG);
%%% The following three options are only available with GASM=1, and in this
%%% order. Options are ignored when inappropriate, see manual. Note that
%%% constants for Sutherland's Law and the cp(T) polynomial are not
%%% available for all fluids
OPT.COMP = 1;   	% 0=incompressible flow; 1=compressible flow (standard); only available with GASM=1 and only valid for gases and incompressible input data; set FLD.C in fluid file
OPT.CPRN = 1;    	% 0=constant Pr-number (input) and k=f(Pr,mu,cp); 1=variable Pr-number: Pr=f(mu(Suth(T)),k(Suth(T)),cp=constant); only available with GASM=1
OPT.CCRP = 1;    	% 0=constant Chapman-Rubesin parameter, 1=general/variable Chapman-Rubesin parameter; only available with GASM=1
%%%
OPT.CPRT = 1;       % 0=constant PrT; 1=variable turbulent Pr-number
OPT.BCIE = 3;       % BC Input Edge [-] (BL Edge input: 1=uE (UE/UI); 2=MaE; 3=psE (PsE/PtI))
OPT.BCEE = 0;       % Wall BC of EE: 0=adiabatic wall; 1=total enthalpy ratio [-]; 2=derivative of total enthalpy ratio [-]; 3=Temperature ratio Twall/TtI [-]; 4=heat flux [W/m2] (heat flux as input only is inaccurate, since unknown wall temperature needs to be estimated)
OPT.SMTH = 0;       % Smoothing experimental data, entry is number of runs)
OPT.TRME = 0;    	% transition method: 0=no transition; 1=prescribed transition location NTR; 2=Wazzan correlation; 3=Michel's Method (Cebeci (2002)); 4=fully turbulent (for air only!)
OPT.RLAM = 0;   	% Relaminarization method, 0=no relam.; 1=prescribed location NRL; 2=simple engineering estimate based on exp. data for air only!
OPT.GRAD = 1;   	% 1=Lagrange; 2=Weighted-difference technique; 3=SPLINE
OPT.CHRT = 0;      	% Pv-diagram and Ts-diagram including range

%% PHYSICAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fluid possibilities: AIR, CO2, MDM2, D6, ...
% cd('../../Fluids')
% MM
% Molecular Fluid properties
FLD.R = 8.3145;                     % [J/mol/K], Universal gas constant
FLD.MW = 162.37752;                  % [g/mol = kg/kmol], Molecular Weight MM ('RefProp')
FLD.Rsg = 1000*FLD.R/FLD.MW;        % [J/kg/K], Specific Gas constant Air
% Fluid Property data (for initial value estimate)
FLD.gamma = 1.0119; %1.011887888336337; % StanMix @T=300 celsius, P=0.01 bar.
FLD.Cp = 1000*FLD.Rsg*FLD.gamma/(FLD.gamma - 1); % [J/kg/K]
FLD.Pr = 0.75;   % [-], Prandlt-number (same IG-limit), StanMix; REFPROP estimate: 0.75;
% FluidProp
FLD.Model = 'StanMix';
FLD.nCmp  = 1;
FLD.Cmp   = 'MM';
FLD.Cnc   = [1 0];
% cd('../INPUT')

%% Case Specific Input Data
load('../Data/ORCHID_Nozzle_Data_StanMix')          % from Excel-file of 10-03-2020 (13-03-2020) 
% load('../Data/ORCHID_PsE_adjusted_StanMix')         % Adjusted PsE (straigthened inlet and outlet, and smoothed kinks at their edges)
load('../Data/ORCHID_PsE_Interpolated_RefProp')  % Use same for both models; interpolated PsE (straigthened inlet and outlet, and smoothed kinks at their edges)

% Geometry (scaled)
INP.L = 1;                                          % [-], scale factor: chord length/plate length/nozzle throat half-width; equal to one for non-scaled geometry
INP.x = transpose(table2array(ORCHID_Nozzle_Data_StanMix(1:249,2))); % [-], from Excel-file surface_flow (10-03-2020)
INP.y = transpose(table2array(ORCHID_Nozzle_Data_StanMix(1:249,3))); % [-], from Excel-file surface_flow (10-03-2020)

%% Free-Stream and BL Edge conditions
INP.PtI = 18.360000e5;                              % [Pa], from Excel-file ORCHID_PR025 (24-01-2020)
INP.TtI = 525.85;                                   % [K], from Excel-file ORCHID_PR025 (24-01-2020)
INP.MaI = table2array(ORCHID_Nozzle_Data_StanMix(1,6));         	% [-], from Excel-file ORCHID_PR025 (24-01-2020)
% INP.PsE = transpose(table2array(ORCHID_Nozzle_Data_StanMix(:,4)));% [-], pE=PsE/PtI, from Excel-file ORCHID_PR025 (24-01-2020)
% INP.PsE = ORCHID_PsE_Interpolated_RefProp; % Use same for both models;
INP.psE = 1/INP.PtI*ORCHID_PsE_Interpolated_RefProp; % [-]

% Wall Boundary Condition (for specific wall enthalpy ratio, derivative of wall enthalpy ratio, wall
% temperature, or wall heat flux (other then adiabatic flow)
if OPT.BCEE > 0
    INP.BCW = zeros(1,length(INP.x));               % [-], scaled h ratio or it's derivative
    % optional
    INP.AWD = [];                                   % only g_aw or T_aw ((static) ratio? -> see OUTPUT) should be given here, calculated after adiabatic simulation
else
    INP.BCW = zeros(1,length(INP.x));               % [-], zero wall heat flux: derivative at wall of enthalpy ratio equal to zero
end

%% Experimental data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%% TURBULENCE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% closure coeficients, for the CS-model (Cebeci, 2002)

% Momentum Equation; and,
TCC.kappa = 0.40;           % [-], Von Karman constant
TCC.A_plus = 26.0;          % [-], Van Driest damping factor constant
TCC.alpha = 0.0168;         % [-], Clauser(’s)/outer eddy viscosity constant
TCC.ints = 11.8;            % [-], Assumed intersection of viscous sublayer with (intermediate) log layer

% Energy Equation:
TCC.kappa_h = 0.44;         % [-], Heat transfer mixing-length constant
TCC.PrT = 0.9;              % [-], (Cebeci (2004)!) Constant Turbulent Prandtl-number (independent of fluid!) if PrT model is not used
TCC.Bcoeff = [34.96 28.79 33.95 6.33 -1.186]; % [-], (Cebeci (1974)!) Turbulent Prandtl-number coefficients

%% NUMERICAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set size of vectors for initilization (vertical/eta-direction)
SET.NPT = 201;                   % preallocation (size is 201 in FORTRAN code)
SET.ITMAX = 6;                   % maximum number of iterations (6 in MAIN code, 8 for IVPL in FORTRAN code)
SET.ITMAX0 = 8; % 20;10;8;       % max number of iterations for initial profile in IVPL (original = 8)
SET.ITMAX2 = 10;                 % maximum number of iterations in PRECAL-file gas model 2
SET.ITMAX3 = 10;                 % maximum number of iterations in PRECAL-file gas model 3
% SET.ITMAX4 = 10;                 % maximum number of iterations in FLDPRS-file (gas model 2)
SET.NTR = []; % 2;               % Force transition at station NTR
SET.NRL = [];                    % Force Re-Laminarization at station NRL

% Grid Parameters
GRD.etaE = 8.0;                  % non-dimensional BL edge location
GRD.VGP = 1.14;                  % variable grid parameter (spacing ratio)
GRD.Deta = 0.01;                 % first grid spacing: h1

%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the following profiles (velocity, dimensionless velocity, enthalpy)
PLT.NSP = [];                    % Plot profiles at stations NS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%