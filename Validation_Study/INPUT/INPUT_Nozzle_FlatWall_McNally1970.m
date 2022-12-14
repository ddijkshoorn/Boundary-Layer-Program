%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT File                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Written by D.D.Dijkshoorn (DDD) on 08/05/2020                         %
%   Checked by ... on ././.                                               %
%   Last checked:  05/03/2022                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Case: NACA0012
%   Source: Fortran program for calculating compressible laminar and turbulent
%   boundary layers in arbitrary pressure gradients, McNally (1970)
%   Validation of: laminar and turbulent flow including point of transition

Case_name = 'Nozzle_FlatWall_McNally1970'; % tag

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
OPT.BCIE = 3;                   % BC Input Edge [-] (BL Edge input: 1=uE (UE/UI); 2=MaE; 3=psE (PsE/PtI))
OPT.BCEE = 3;0;                   % Wall BC of EE: 0=adiabatic wall; 1=total enthalpy ratio [-]; 2=derivative of total enthalpy ratio [-]; 3=Temperature ratio Twall/TtI [-]; 4=heat flux [W/m2] (heat flux as input only is inaccurate, since unknown wall temperature needs to be estimated)
OPT.SMTH = 0;                   % Smoothing experimental data, entry is number of runs)
OPT.TRME = 4;                   % transition method: 0=no transition; 1=prescribed transition location NTR; 2=Wazzan correlation; 3=Michel's Method (Cebeci (2002)); 4=fully turbulent
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
INP.L = 1;                      % [m], scale factor: chord length/plate length/nozzle throat half-width; equal to one for non-scaled geometry
INP.x = ...                     % [m]
0.3048*[0       0.08330 0.16700 0.25000 0.33300 0.41670 0.50000 0.58330 0.66700...
        0.75000 0.83300 0.91670 1.00000 1.08330 1.16700];

INP.y = zeros(1,length(INP.x)); % [m]

%% Free-Stream and BL Edge conditions
INP.PtI = 47.880259*2004.77;        % [Pa], pound per square foot to Pascal
INP.TtI = 532.00/1.8000;            % [K], degrees Rankine to Kelvin
INP.MaI = 0.4821;               % [-]
INP.delta_astI = 0.3048*0.004710;% [m]
INP.thetaI = 0.3048*0.003420;   % [m]
% INP.PsE = ... % [Pa], obtained with grabit
% INP.PtI*[0.850023464658200 0.849521173208322 0.845668314732295 0.841386331955570...
%          0.834992649478984 0.826860393243863 0.815491098278653 0.801429320909659...
%          0.781411158647631 0.751097660459126 0.706992361650832 0.653709559010792...
%          0.594903323442442 0.532844187504718 0.466159304693922];
INP.psE = ... % [Pa], obtained with grabit (changed to psE on 05/03/2022)
        [0.850023464658200 0.849521173208322 0.845668314732295 0.841386331955570...
         0.834992649478984 0.826860393243863 0.815491098278653 0.801429320909659...
         0.781411158647631 0.751097660459126 0.706992361650832 0.653709559010792...
         0.594903323442442 0.532844187504718 0.466159304693922];

% Wall Boundary Condition (for specific wall enthalpy ratio, derivative of wall enthalpy ratio, wall
% temperature, or wall heat flux (other then adiabatic flow)
if OPT.BCEE > 0
    INP.BCW = ones(1,length(INP.x));   % [-], Tw/TtI = 1
    INP.AWD = []; % T_aw [K]
else
    INP.BCW = zeros(1,length(INP.x));   % [-], zero wall heat flux: derivative at wall of enthalpy ratio equal to zero
end

% % Interpolation
% xq = [0:0.001:INP.x(end) INP.x(end)];
% p = pchip(INP.x,INP.y,xq);
% q = pchip(INP.x,INP.PsE,xq);
% INP.x = xq;
% INP.y = p;
% INP.PsE = q;
% %

%% Corrected for initial BL thickness (by adding a flat plate of length L_FP
%%% L_FP = 0.6850 m, size determined with other calculation file by avering
%%% x-coordinate between correct delta_ast and theta value.
INP.L_FP = 0.6850;
INP.dx0 = 0.005;
INP.x0 = 0:INP.dx0:INP.L_FP;
INP.x = [INP.x0 (INP.x(2:end) + INP.L_FP)];
INP.y = [INP.y(1)*ones(1,length(INP.x0)) INP.y(2:end)];
% INP.PsE = [INP.PsE(1)*ones(1,length(INP.x0)) INP.PsE(2:end)];
INP.psE = [INP.psE(1)*ones(1,length(INP.x0)) INP.psE(2:end)]; % Changed to pse on 05/03/2022

% INP.BCW = zeros(1,length(INP.x)); % adiabatic wall Temp. calculation
INP.BCW = ones(1,length(INP.x));   % [-], Tw/TtI = 1
% Taw was calculated in separate run
INP.AWD = [293.330818885335,293.738195252866,294.118713827942,294.208781295537,294.243413762481,294.307864445948,294.322013862311,294.354567452795,294.372353925175,294.382728249878,294.401353521494,294.411762134235,294.423376697002,294.431677045654,294.439625750438,294.447781329230,294.454714071029,294.461290902131,294.468773237009,294.470680042896,294.479712817668,294.483556387502,294.490023399739,294.493750108261,294.496807312233,294.502333816679,294.505780248217,294.509285785025,294.512377611083,294.515975527864,294.518572874065,294.521919590884,294.524562231467,294.528288644272,294.530522157762,294.532370217875,294.535607460442,294.537800364982,294.540344540631,294.542163859081,294.544422898339,294.546297736373,294.548207837694,294.550042133701,294.551886890581,294.553584106294,294.555416009992,294.557061194894,294.558809918387,294.560955886606,294.562454261367,294.563574768651,294.565567407615,294.566913231495,294.568686977108,294.569916017972,294.571393611929,294.572648517231,294.573964740132,294.575127690762,294.576378931212,294.577481157012,294.578651966948,294.579742618820,294.580847980174,294.581917880435,294.582996984996,294.584041365475,294.585101599748,294.586135681642,294.587177032880,294.588230143136,294.589302938699,294.590685670705,294.591623404664,294.592335501772,294.593566699498,294.594393427150,294.595534965055,294.596405894839,294.597375402161,294.598219203057,294.599132407848,294.599923684944,294.600759657702,294.601533053688,294.602304762199,294.603041941713,294.603782974758,294.604484800164,294.605195591005,294.605881541939,294.606560531170,294.607234727024,294.607894761157,294.608554034591,294.609203841645,294.609851084389,294.610491903372,294.611133711103,294.611769962515,294.612426037994,294.613084047403,294.613758739427,294.614411112470,294.615070686251,294.615727985722,294.616575625202,294.617160046312,294.617607385228,294.618384275493,294.618893882580,294.619604134272,294.620177093091,294.620824256746,294.621381648035,294.621988595509,294.622539803115,294.623104610026,294.623632753126,294.624172545663,294.624676929197,294.625186876694,294.625676904815,294.626159252524,294.626631339728,294.627095927711,294.627549771561,294.627998211383,294.628439629573,294.628871937426,294.629303877162,294.629723871300,294.630145949702,294.630557980395,294.630970905149,294.631376530950,294.631667464859,294.618764798050,294.578794732324,294.530301599850,294.464632772788,294.378554186203,294.262718430618,294.114125161788,293.893371968362,293.549095764113,293.055587756094,292.450249578662,291.775749379257,291.022433909893,290.103237716363];
% INP.AWD = []; %%% added for interpolation

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
SET.NTR = [];                    % Force transition at station NTR

% Grid Parameters
GRD.etaE = 8.0;                  % non-dimensional BL edge location
GRD.VGP = 1.14; % 1;1.02; for smooth BL velocity thickness                 % variable grid parameter (spacing ratio)
GRD.Deta = 0.01;     % 0.1;      % first grid spacing: h1

%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the following profiles (velocity, dimensionless velocity, enthalpy)
PLT.NSP = [];                   % Plot profiles at stations NS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%