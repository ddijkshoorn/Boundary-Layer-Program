%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT File                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Written by D.D.Dijkshoorn (DDD) on 13/06/2019             %
%   Checked by ... on ././.                                   %
%   Last revision: ...                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Case: Flat Plate
%   Source: to be determined
%   Verification of: to be determind

Case_name = 'FlatPlate'; % tag

%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ideal Gas (IG); Real Gas (RG)

OPT.CASE = 1;                   % Case number: 1=NACA0012 Cebeci(2002) (for results output in PLOTFILE and TABLEGEN)
OPT.GASM = 1;                   % Gas model: 1=calorically perfect IG; 2=thermally perfect IG; 3=Real Gas (RG)
OPT.COMP = 1;                   % Compressible flow (standard), incompressible=0 only available with GASM=1
OPT.CCRP = 0;    %%% !!! %%%    % 0=constant Chapman-Rubesin parameter, 1=general/variable Chapman-Rubesin parameter
OPT.CPRN = 0;
OPT.CPRT = 0;                   % 0=constant PrT; 1=variable turbulent Pr-number
OPT.BCIE = 2;                   % BC Input Edge (BL Edge input: 1=UE/uE; 2=MaE; 3=psE (ratio?))
OPT.SMTH = 0;                   % Smoothing experimental data, entry is number of runs)
OPT.GRAD = 1;                   % 1=Lagrange; 2=Weighted-difference technique; 3=SPLINE
OPT.BCEE = 1;                   % Wall BC of EE: 0=adiabatic wall; 1=enthalpy ratio [-]; 2=derivative of enthalpy ratio [-]; 3=Temperature ratio [-]; 4=heat flux [W/m2] (temperature or heat flux directly)
OPT.IBLT = 0;                   % 1=Start calculation with initial BL Thickness
OPT.TRME = 0;                   % transition method: 0=no transition; 1=prescribed transition location NTR; 2=Wazzan correlation; 3=Michel's Method (Cebeci (2002)); 4=fully turbulent
OPT.CHRT = 0;                   % Pv-diagram and Ts-diagram including range
OPT.MURS = 0;                   % multiple runs
OPT.SOLV = 1;                   % 1=algebraic solver; 2=matrix solver

%% PHYSICAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Case Specific Input Data

% Geometry (scaled)
INP.L = 1;                      % [m], scale factor: chord length/plate length/nozzle throat half-width; equal to one for non-scaled geometry
INP.x = 0;
INP.y = zeros(1,length(INP.x)); % [-]

%% Fluid Properties

% Fluid (instead of OPT.FLUID) Possibilities: AIR, CO2, MDM2, D6
cd('../../Fluids')
AIR_Rogers1992
cd('..')

% Case specific fluid properties (overriding AIR-file)
FLD.Rsg = 287;                    % [J/kg/K], for air
FLD.gamma = 1.4;                % [-], @300 Celsius (White, 2006)
FLD.Cp = FLD.Rsg*FLD.gamma/(FLD.gamma - 1); % [J/kg/K]
FLD.mu = 2.93e-5;               % [kg/m/s], @300 Celsius (White, 2006)
FLD.k = 4.39e-2;                % [W/m/K], @300 Celsius (White, 2006)
FLD.Pr = 1;0.723;               % [-], @300 Celsius (White, 2006)
FLD.C = 1.0;                    % [-], Predfined Chapman-Rubesin parameter (constant)
% if INP.Pr == 1 %%% Why?
%     INP.MaI = 0;
% end
INP.Tref = 0;

%% Free-Stream and BL Edge conditions
INP.PtI = 7.8244e+05;           % [Pa], for 1 bar static pressure
INP.TtI = 500;                  % [K], for approximately 20 degrees Celsius
INP.MaI = 0.0;                  % [m/s]
INP.MaE = INP.MaI*ones(1,length(INP.x));% [-]

% Wall Boundary Condition (for specific wall enthalpy ratio, derivative of wall enthalpy ratio, wall
% temperature, or wall heat flux (other then adiabatic flow)
% if OPT.BCEE > 0
%     INP.BCW = INP.gw*ones(1,length(INP.x));   % [-], scaled h ratio or it's derivative
%     % optional %%% Check!!! %%%
%     INP.AWD = []; % only g_aw or T_aw ((static) ratio? -> see OUTPUT) should be given here, calculated after adiabatic simulation
% else
%     INP.BCW = zeros(1,length(INP.x));   % [-], zero wall heat flux: derivative at wall of enthalpy ratio equal to zero
% end

%% NUMERICAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set size of vectors for initilization (vertical/eta-direction)
SET.NPT = 1000;201;  % 600;            % preallocation (size is 201 in FORTRAN code)
SET.ITMAX = 20;10;6;                   % maximum number of iterations (6 in MAIN code, 8 for IVPL in FORTRAN code)
SET.NTR = [];                    % Force transition at station NTR

% Grid Parameters - adapted for laminar flows to constant VGP
GRD.etaE = 8.0;                  % non-dimensional BL edge location
GRD.VGP = 1.0;1.14; % 1;1.02; for smooth BL velocity thickness                 % variable grid parameter (spacing ratio)
GRD.Deta = 0.01;0.01;     % 0.1;      % first grid spacing: h1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% used for accurate results at the wall (f''(0)) with m=1000 and/or gw=0.0:
SET.NPT = 100000;
% SET.ITMAX = 6;
GRD.etaE = 8.0;
GRD.VGP = 1.0;
GRD.Deta = 0.0001;