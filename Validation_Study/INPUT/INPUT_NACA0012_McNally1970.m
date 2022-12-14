%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT File                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Written by D.D.Dijkshoorn (DDD) on 08/05/2020                         %
%   Checked by ... on ././.                                               %
%   Last check: 05/03/2022                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Case: NACA0012
%   Source: Fortran program for calculating compressible laminar and turbulent
%   boundary layers in arbitrary pressure gradients, McNally (1970)
%   Validation of: laminar and turbulent flow including point of transition

Case_name = 'NACA0012_McNally1970'; % tag

%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ideal Gas (IG); Nonideal Gas (NG)

OPT.GASM = 1;                   % Gas model: 1=calorically perfect IG; 2=thermally perfect IG; 3=Nonideal Gas (NG);
%%% The following three options are only available with GASM=1, and in this
%%% order. Options are ignored when inappropriate, see manual. Note that
%%% constants for Sutherland's Law and the cp(T) polynomial are not
%%% available for all fluids
OPT.COMP = 1;                   % 0=incompressible flow; 1=compressible flow (standard); only available with GASM=1 and only valid for gases and incompressible input data; set FLD.C in fluid file
OPT.CPRN = 0;                   % 0=constant Pr-number (input) and k=f(Pr,mu,cp); 1=variable Pr-number: Pr=f(mu(Suth(T)),k(Suth(T)),cp=constant); only available with GASM=1
OPT.CCRP = 1;                   % 0=constant Chapman-Rubesin parameter, 1=general/variable Chapman-Rubesin parameter; only available with GASM=1
%%%
OPT.CPRT = 0;                   % 0=constant PrT; 1=variable turbulent Pr-number
OPT.BCIE = 1;                   % BC Input Edge (BL Edge input: 1=UE/uE; 2=MaE; 3=psE (ratio?))
OPT.BCEE = 3;0;                   % Wall BC of EE: 0=adiabatic wall; 1=enthalpy ratio [-]; 2=derivative of enthalpy ratio [-]; 3=Temperature ratio [-]; 4=heat flux [W/m2] (temperature or heat flux directly)
OPT.SMTH = 0;                   % Smoothing experimental data, entry is number of runs)
OPT.TRME = 1;3;                   % transition method: 0=no transition; 1=prescribed transition location NTR; 2=Wazzan correlation; 3=Michel's Method (Cebeci (2002)); 4=fully turbulent
OPT.RLAM = 0;                   % 1=Relaminarization check, 0=no check; simple engineering estimate based on exp. data for air
OPT.GRAD = 1;                   % 1=Lagrange; 2=Weighted-difference technique; 3=SPLINE
OPT.CHRT = 0;                   % Pv-diagram and Ts-diagram including range

%% PHYSICAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fluid possibilities: AIR, CO2, MDM2, D6, ...
% cd('../Fluids')
% % AIR
% AIR_NACA0012_Cebeci2002
% cd('../INPUT')

% Molecular Fluid properties
FLD.R = 8.3145;                 % [J/mol/K], Universal gas constant
FLD.MW = 28.97;                 % [g/mol = kg/kmol], Molecular Weight of dry Air
FLD.Rsg = FLD.R/FLD.MW*1e3;     % [J/kg/K], Specific Gas constant Air

% Fluid Property data
FLD.gamma = 1.4;                % [-], Specific heat ratio (Cp/Cv)
FLD.SLV = 110.4;                % [K],  (Groot, 2018) Sutherland's constant (later ones: White (2006) and Cebeci (2002)!)
FLD.SLC = 194;                  % [K], (Groot, 2018; and White, 2006) Sutherland's constant for thermal conductivity
FLD.Tref_S = 273;               % [K], (Groot, 2018; and White, 2006) reference temperature Sutherland's Law
FLD.mu_ref = 1.716e-5;          % [Pas], (Groot, 2018; and White, 2006) reference viscosity Sutherland's Law
FLD.k_ref = 0.0241;             % [W/m/K], (Groot, 2018; and White, 2006) reference thermal conductivity Sutherland's Law
FLD.Pr = 0.69750;0.72;          % [-], Prandlt-number
FLD.Cp = FLD.Rsg*FLD.gamma/(FLD.gamma - 1); % [J/kg/K]
% FLD.Cpcoeff = [1022.5294853 -0.1758625 4.020605e-4 -4.8640623e-8]; % a0, a1, a2, a3, a4, etc. -> Cp = a0 + a1*T + ... % Andrews (1981) IG relation; range: 100 - 590 K

%% Case Specific Input Data

% Geometry (scaled)
INP.L = 1;   %%% are you sure? CHECK!!!      % [m], scale factor: chord length/plate length/nozzle throat half-width; equal to one for non-scaled geometry
INP.x = ...                     % [-]
0.3048*[0.      0.02500 0.06250 0.12500 0.25000 0.37500 0.50000 0.75000...
        1.00000 1.25000 1.30000 1.35000 1.40000 1.45000 1.50000 1.55000...
        1.60000 1.65000 1.70000 1.75000 2.00000 2.25000 2.50000 2.75000...
        3.00000 3.25000 3.50000 3.75000 4.00000 4.25000 4.50000 4.75000 5.00000];
INP.y = ...
0.3048*[0.      0.05765 0.09470 0.13075 0.17775 0.21000 0.23415 0.26725...
        0.28685 0.29705 0.29810 0.29890 0.29960 0.29995 0.30010 0.29995...
        0.29960 0.29900 0.29825 0.29735 0.29015 0.27905 0.26470 0.24760...
        0.22815 0.20685 0.18320 0.15805 0.13115 0.10275 0.07240 0.04035...
        0.00630];

%% Free-Stream and BL Edge conditions
INP.PtI = 47.880259*2645.00;    % [Pa], pound per square feet to Pascal
INP.TtI = 600/1.800;            % [K], degrees Rankine to Kelvin
INP.MaI = 0.2840;               % [m/s], feet per second to meter per second
%%% calculated (here instead of PRECAL-file because of non-matching edge data)
INP.TsI = INP.TtI/(1 + (FLD.gamma - 1)/2*INP.MaI^2);            % [K]
aI = sqrt(FLD.gamma*FLD.Rsg*INP.TsI);                           % [m/s]
INP.UI = INP.MaI*aI;                                            % [m/s]
INP.PsI = INP.PtI*(INP.TsI/INP.TtI)^(FLD.gamma/(FLD.gamma - 1));% [Pa]
%%%
INP.uE = ...                    % [-]
0.3048/INP.UI*[0.      270.64000 339.99200 376.86600 397.16400 400.54700 401.90000...
      401.90000 400.20900 397.16400 396.38600 395.57400 394.79600 393.95000...
      393.10500 392.22500 391.31200 390.39800 389.48500 388.53800 383.97100...
      379.43700 374.83600 370.10000 365.36400 360.89800 356.23000 351.15500...
      345.74300 338.63800 330.85700 322.06200 311.91300];

% Wall Boundary Condition (for specific wall enthalpy ratio, derivative of wall enthalpy ratio, wall
% temperature, or wall heat flux (other then adiabatic flow)
if OPT.BCEE > 0
    INP.BCW = ones(1,length(INP.x));   % [-], Tw/TtI = 1
    % Taw calculated with separate run
    INP.AWD = [333.333333333333,332.746549288760,332.427748478791,332.229182412749,332.108010447885,332.091128761979,332.088305406848,332.093240280592,332.105157274726,332.126184825275,332.143308498202,332.199936812984,332.304858692667,332.429107585712,332.539885337735,332.624413080625,332.676482225264,332.714255717999,332.745693888303,332.761866142900,332.814969336697,332.846859113516,332.855543598299,332.891952652945,332.900768826030,332.932794701473,332.934247097382,332.965534939684,332.971869041021,333.005373631159,333.018123628364,333.053051402046,333.070372613550];
else
    INP.BCW = zeros(1,length(INP.x));   % [-], zero wall heat flux: derivative at wall of enthalpy ratio equal to zero
end

%% Experimental data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% included in specific PLOT-file for generating thesis graphs
% cd('./Data')
% load('NACA0012_FORTRAN_SimResults_Cebeci2002')
% cd('..')

%% TURBULENCE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Closure coeficients, for the CS-model: EV (Cebeci, 2002), and PrT (Cebeci, 1974)

% Momentum Equation; and,
TCC.kappa = 0.40;           % [-], Von Karman constant
TCC.A_plus = 26.0;          % [-], Van Driest damping factor constant
TCC.alpha = 0.0168;         % [-], Clauser(’s)/outer eddy viscosity constant
TCC.ints = 11.8;            % [-], Assumed intersection of viscous sublayer with (intermediate) log layer

% Energy Equation:
TCC.kappa_h = 0.44;         % [-], Heat transfer mixing-length constant
TCC.PrT = 0.9;              % [-], (Cebeci (2004)!) Constant Turbulent Prandtl-number (independent of fluid!) if PrT model is not used
TCC.Bcoeff = [34.96 28.79 33.95 6.33 -1.186]; % [-], (Cebeci (1974)!) Turbulent Prandtl-number coefficients

%% NUMERICAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set size of vectors for initilization (vertical/eta-direction)
SET.NPT = 201;  % 600;           % preallocation (size is 201 in FORTRAN code)
SET.ITMAX = 6;                   % maximum number of iterations (6 in MAIN code, 8 for IVPL in FORTRAN code)
SET.ITMAX0 = 8; % 20;10;8;       % max number of iterations for initial profile in IVPL (original = 8)
% SET.ITMAX3 = 10;                 % maximum number of iterations (PRECAL-file gas model 2)
SET.NTR = 10;8; % 10 = perfect!!   % Force transition at station NTR

% Grid Parameters
GRD.etaE = 8.0;                  % non-dimensional BL edge location
GRD.VGP = 1.14; % 1;1.02; for smooth BL velocity thickness                 % variable grid parameter (spacing ratio)
GRD.Deta = 0.01;     % 0.1;      % first grid spacing: h1

%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the following profiles (velocity, dimensionless velocity, enthalpy)
PLT.NSP = [];                   % Plot profiles at stations NS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%