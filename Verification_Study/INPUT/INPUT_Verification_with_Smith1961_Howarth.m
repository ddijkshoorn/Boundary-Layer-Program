%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT File                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Written by D.D.Dijkshoorn (DDD) on 19/04/2018             %
%   Checked by ... on ././.                                   %
%   Last revision: ...                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Case: NACA0012
%   Source: Convective Heat Transfer, Cebeci (2002)
%   Verification of: solver, BL variables/characteristics,
%   turbulence model

Case_name = 'HowarthsFlow_Smith1961'; % tag

%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ideal Gas (IG); Real Gas (RG)

OPT.CASE = 1;                   % Case number: 1=NACA0012 Cebeci(2002) (for results output in PLOTFILE and TABLEGEN)
OPT.GASM = 1;                   % Gas model: 1=calorically perfect IG; 2=thermally perfect IG; 3=Real Gas (RG)
OPT.COMP = 0;                   % Compressible flow (standard), incompressible=0 only available with GASM=1
OPT.CCRP = 1;    %%% !!! %%%    % 0=constant Chapman-Rubesin parameter, 1=general/variable Chapman-Rubesin parameter
OPT.CPRN = 0;
OPT.CPRT = 0;                   % 0=constant PrT; 1=variable turbulent Pr-number
OPT.BCIE = 1;                   % BC Input Edge (BL Edge input: 1=UE/uE; 2=MaE; 3=psE (ratio?))
OPT.SMTH = 0;                   % Smoothing experimental data, entry is number of runs)
OPT.GRAD = 1;                   % 1=Lagrange; 2=SPLINE; 3=Weighted-difference technique
OPT.BCEE = 0;                   % Wall BC of EE: 0=adiabatic wall; 1=enthalpy ratio [-]; 2=derivative of enthalpy ratio [-]; 3=Temperature ratio [-]; 4=heat flux [W/m2] (temperature or heat flux directly)
OPT.IBLT = 0;                   % 1=Start calculation with initial BL Thickness
OPT.TRME = 0;                   % transition method: 0=no transition; 1=prescribed transition location NTR; 2=Wazzan correlation; 3=Michel's Method (Cebeci (2002)); 4=fully turbulent
OPT.CHRT = 0;                   % Pv-diagram and Ts-diagram including range
OPT.MURS = 0;                   % multiple runs
OPT.SOLV = 1;                   % 1=algebraic solver; 2=matrix solver

%% PHYSICAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fluid (instead of OPT.FLUID) Possibilities: AIR, CO2, MDM2, D6
cd('../../Fluids')
AIR
cd('../Verification_Study')
% cd('..')

%% Case Specific Input Data

% Geometry (scaled)
INP.L = 1;   %%% are you sure? CHECK!!!      % [m], scale factor: chord length/plate length/nozzle throat half-width; equal to one for non-scaled geometry
% NB I added a zero at X(0):
INP.x = [0 0.0125 0.025 0.050 0.075 0.100 0.150 0.200 0.300 0.400 0.600 0.800 0.840 0.880 0.920 0.948 0.956 0.958 0.9589];                    % [-]
% INP.x = 0:0.0001:1.0; % [-]
INP.y = zeros(1,length(INP.x));% [-]

%% Free-Stream and BL Edge conditions
INP.PsI = 1e5;                  % [Pa]
INP.TsI = 500;                  % [K]
INP.UI = 1 - 1/8*INP.x(1);      % [m/s]
INP.uE = (1 - 1/8*INP.x)./INP.UI;  % [m/s], uE = UE/U0?
% niet ideale oplossing hierboven

% Wall Boundary Condition (for specific wall enthalpy ratio, derivative of wall enthalpy ratio, wall
% temperature, or wall heat flux (other then adiabatic flow)
if OPT.BCEE > 0
    INP.BCW = zeros(1,length(INP.x));   % [-], scaled h ratio or it's derivative
    % optional %%% Check!!! %%%
    INP.AWD = []; % only g_aw or T_aw ((static) ratio? -> see OUTPUT) should be given here, calculated after adiabatic simulation
else
    INP.BCW = zeros(1,length(INP.x));   % [-], zero wall heat flux: derivative at wall of enthalpy ratio equal to zero
end

%% TURBULENCE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% closure coeficients, for the

% Momentum Equation; and,
TCC.kappa = 0.40;           % [-], Von Karman constant
TCC.A_plus = 26.0;          % [-], Van Driest damping factor constant
TCC.alpha = 0.0168;         % [-], Clauser(’s)/outer eddy viscosity constant
TCC.ints = 11.8;            % [-], Assumed intersection of viscous sublayer with (intermediate) log layer

% Energy Equation:
TCC.kappa_h = 0.44;         % [-], Heat transfer mixing-length constant
TCC.PrT = 0.9;              % [-], (Cebeci (2004)!) Constant Turbulent Prandtl-number (independent of fluid!) if PrT model is not used
TCC.Bcoeff = [34.96 28.79 33.95 6.33 -1.186]; % [-], (Cebeci (1974)!) Turbulent Prandtl-number coefficients

%% NUMERICAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set size of vectors for initilization (vertical/eta-direction)
SET.NPT = 1000;%100000;              % preallocation (size is 201 in FORTRAN code)
SET.ITMAX = 6;              % maximum number of iterations (6 in MAIN code, 8 for IVPL in FORTRAN code)
SET.NTR = [];               % Force transition at station NTR

% Grid Parameters
GRD.etaE = 9.0;             % non-dimensional BL edge location
GRD.VGP = 1.0;              % variable grid parameter (spacing ratio)
GRD.Deta = 0.01;%0.0001;            % first grid spacing: h1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%