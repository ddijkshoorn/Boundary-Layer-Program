%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT File                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Written by D.D.Dijkshoorn (DDD) on 23/09/2017                         %
%   Checked by ... on ././.                                               %
%   Last revision: 26/11/2022                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Case: 
%   Source: 
%   Verification/Validation of: 

Case_name = 'XXX'; % tag, NACA0012 (Cebeci, 2002) taken as example here

%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ideal Gas (IG); Nonideal Gas (NG)

OPT.GASM = 1;                   % Gas model: 1=calorically perfect IG; 2=thermally perfect IG; 3=Nonideal Gas (NG);
%%% The following three options are only available with GASM=1, and in this
%%% order. Options are ignored when inappropriate, see manual. Note that
%%% constants for Sutherland's Law and the cp(T) polynomial are not
%%% available for all fluids
OPT.COMP = 1;                   % Compressible flow (standard), incompressible=0; only available with GASM=1 and only valid for gases and incompressible input data
OPT.CPRN = 0;                   % 0=constant Pr-number (input) and k=f(Pr,mu,cp); 1=variable Pr-number: Pr=f(mu(Suth(T)),k(Suth(T)),cp); only available with GASM=1
OPT.CCRP = 1;                   % 0=constant Chapman-Rubesin parameter, 1=general/variable Chapman-Rubesin parameter; only available with GASM=1
%%%
OPT.CPRT = 0;                   % 0=constant PrT; 1=variable turbulent Pr-number
OPT.BCIE = 1;                   % BC Input Edge [-] (BL Edge input: 1=uE (UE/UI); 2=MaE; 3=psE (PsE/PtI))
OPT.BCEE = 0;                   % Wall BC of EE: 0=adiabatic wall; 1=total enthalpy ratio [-]; 2=derivative of total enthalpy ratio [-]; 3=Temperature ratio Twall/TtI [-]; 4=heat flux [W/m2] (heat flux as input only is inaccurate, since unknown wall temperature needs to be estimated)
OPT.SMTH = 0;                   % Smoothing experimental data, entry is number of runs)
OPT.TRME = 1;                   % Transition method: 0=no transition; 1=prescribed transition location NTR; 2=Wazzan correlation; 3=Michel's Method (Cebeci (2002)); 4=fully turbulent (for air only!)
OPT.RLAM = 0;                   % Relaminarization method, 0=no relam.; 1=prescribed location NRL; 2=simple engineering estimate based on exp. data for air only!
OPT.GRAD = 1;                   % 1=Lagrange; 2=Weighted-difference technique; 3=SPLINE
OPT.CHRT = 0;                   % Pv-diagram and Ts-diagram including range (not available)

%% PHYSICAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fluid possibilities: AIR, CO2, MDM2, D6, ...
cd('../Fluids')
AIR
cd('../INPUT')

%% Case Specific Input Data

% Geometry (scaled)
INP.L = 1;                      % [m], scale factor: chord length/plate length/nozzle throat half-width; equal to one for non-scaled geometry
INP.x = [];                     % [-]
INP.y = [];                     % [-]

%% Free-Stream and BL Edge conditions (3 options) (remove rest, or leave unused option empty)
% 1) Experimental data (smooth if needed): velocity distribution
if OPT.BCIE == 1 % can be removed, example
INP.PsI = 1.2e5;                % [Pa]
INP.TsI = 300;                  % [K]
INP.UI = [];                    % [m/s]
INP.uE = ones(1,33);            % [-]

% 2) Experimental simulations: Mach distribution
elseif OPT.BCIE == 2 % can be removed
INP.PtI = [];                   % [Pa]
INP.TtI = [];                   % [K]
INP.MaI = [];                   % [-]
INP.MaE = [];                   % [-]

% 3) Nozzle/stator (scaled): pressure distribution
elseif OPT.BCIE == 3 % can be removed
INP.PtI = [];                   % [Pa]
INP.TtI = [];                   % [K]
INP.MaI = [];                   % [-]
INP.psE = [];                   % [-], psE=PsE/PtI: static (surface) pressure scaled with PtI
end

% Wall Boundary Condition (for specific wall enthalpy ratio, derivative of wall enthalpy ratio, wall
% temperature, or wall heat flux (other then adiabatic flow)
if OPT.BCEE > 0
    INP.BCW = zeros(1,length(INP.x));   % [-], scaled h ratio or it's derivative
    % optional
    INP.AWD = []; % only g_aw (for OPT.BCEE = 1 or 2) or T_aw [K] (for OPT.BCEE = 3 or 4), which can be calculated with adiabatic simulation first
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

%% NUMERICAL Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set size of vectors for initilization (vertical/eta-direction)
SET.NPT = 201;  % 600;           % preallocation (size is 201 in original FORTRAN code)
SET.ITMAX = 6;                   % maximum number of iterations (6 in MAIN code, 8 for IVPL in FORTRAN code)
SET.ITMAX0 = 8; % 20;10;8;       % max number of iterations for initial profile in IVPL (original = 8)
SET.ITMAX2 = 10;                 % maximum number of iterations in PRECAL-file gas model 2
SET.ITMAX3 = 10;                 % maximum number of iterations in PRECAL-file gas model 3
% SET.ITMAX4 = 10;                 % maximum number of iterations in FLDPRS-file (gas model 2)
SET.NTR = 53;                    % Force transition at station NTR
SET.NRL = [];                    % Force Re-Laminarization at station NRL

% Grid Parameters
GRD.etaE = 8.0;                  % non-dimensional BL edge location
GRD.VGP = 1.14; % 1;1.02; for smooth BL velocity thickness    % variable grid parameter (spacing ratio)
GRD.Deta = 0.01;     % 0.1;      % first grid spacing: h1

%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the following profiles (velocity, dimensionless velocity, enthalpy)
PLT.NSP = [23 55 73 87];               % Plot profiles at stations NS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%