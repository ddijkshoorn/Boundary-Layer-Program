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

Case_name = 'Nozzle_CurvedWall_McNally1970'; % tag

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
OPT.RLAM = 0;1;                   % 1=Relaminarization check, 0=no check; 1=simple engineering estimate based on exp. data for air
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
0.3048*[0       0.08330 0.16700 0.25000 0.33300 0.41670 0.50000 0.58330...
        0.66700 0.75000 0.83300 0.91670 1.00000 1.08330 1.16700];
INP.y = ...
0.3048*[0.58300 0.58300 0.58300 0.58300 0.58300 0.58300 0.58300 0.57300...
        0.54300 0.49700 0.45100 0.42500 0.41670 0.42500 0.45100];

%% Free-Stream and BL Edge conditions
INP.PtI = 47.880259*2004.77;        % [Pa], pound per square foot to Pascal
INP.TtI = 532.00/1.8000;            % [K], degrees Rankine to Kelvin
INP.MaI = 0.4821;               % [-]
INP.delta_astI = 0.3048*0.004710;% [m]
INP.thetaI = 0.3048*0.003420;   % [m]
% INP.PsE = ... % [Pa]
% INP.PtI*[0.853000 0.854300 0.855100 0.858300 0.862300 0.871400 0.889600...
%          0.896800 0.882000 0.800000 0.653000 0.498000 0.355000 0.235000...
%          0.124000];
INP.psE = ... % [Pa]
        [0.853000 0.854300 0.855100 0.858300 0.862300 0.871400 0.889600...
         0.896800 0.882000 0.800000 0.653000 0.498000 0.355000 0.235000...
         0.124000]; % changed to pse on 5/03/2022

% % Interpolation
% xq = [0:0.001:INP.x(end) INP.x(end)];
% p = pchip(INP.x,INP.y,xq);
% q = pchip(INP.x,INP.PsE,xq);
% INP.x = xq;
% INP.y = p;
% INP.PsE = q;
% %

% Wall Boundary Condition (for specific wall enthalpy ratio, derivative of wall enthalpy ratio, wall
% temperature, or wall heat flux (other then adiabatic flow)
if OPT.BCEE > 0
    INP.BCW = ones(1,length(INP.x));   % [-], Tw/TtI = 1
    % Taw was calculated in separate run
    INP.AWD = [293.377743098174,293.980652468847,294.453808576291,294.457663278655,294.578938764824,294.718504943106,294.884168437975,294.913898179317,294.622629568082,293.632598275286,291.917821093754,290.072381012714,287.790047445876,285.212127033870,281.059230451827]; % T_aw [K]
else
    INP.BCW = zeros(1,length(INP.x));   % [-], zero wall heat flux: derivative at wall of enthalpy ratio equal to zero
end

%% Corrected for initial BL thickness (by adding a flat plate of length L_FP
%%% L_FP = 0.6850 m, size determined with other calculation file by avering
%%% x-coordinate between correct delta_ast and theta value.
INP.L_FP = 0.6850;
INP.dx0 = 0.005;
INP.x0 = 0:INP.dx0:INP.L_FP;
INP.x = [INP.x0 (INP.x(2:end) + INP.L_FP)];
INP.y = [INP.y(1)*ones(1,length(INP.x0)) INP.y(2:end)];
% INP.PsE = [INP.PsE(1)*ones(1,length(INP.x0)) INP.PsE(2:end)];
INP.psE = [INP.psE(1)*ones(1,length(INP.x0)) INP.psE(2:end)]; % changed to pse on 05/03/2022

% % % INP.BCW = zeros(1,length(INP.x)); % adiabatic wall Temp. calculation
INP.BCW = ones(1,length(INP.x));   % [-], Tw/TtI = 1
% determined with separate calculation
INP.AWD = [293.377743098174,293.785308330404,294.141870483424,294.230487392163,294.272110167999,294.332011823634,294.347653590564,294.378132231161,294.396491261429,294.406295263858,294.424582356848,294.435100587948,294.446172593840,294.454744647724,294.462096674239,294.470522565587,294.476931725720,294.483765504163,294.490804579437,294.492993586701,294.501550052609,294.505540676753,294.511804108186,294.515384586103,294.518475662773,294.523811312065,294.527283820277,294.530653331776,294.533759958643,294.537221158460,294.539848300123,294.543053595680,294.545717980518,294.549312564793,294.551513999034,294.553291226318,294.556498904923,294.558651091271,294.561135525828,294.562937541047,294.565152449979,294.566992032055,294.568873026283,294.570672390628,294.572482540585,294.574152516037,294.575944851844,294.577565132908,294.579272996597,294.580923669221,294.583049864306,294.584529403234,294.585621827409,294.587540753031,294.588856462986,294.590569522353,294.591786092320,294.593203863577,294.594436610131,294.595711266655,294.596840238350,294.598053283161,294.599126339568,294.600249428238,294.601316233142,294.602372562199,294.603412602761,294.604449310678,294.605456024738,294.606477122831,294.607473959094,294.608481236116,294.609514544072,294.610836148703,294.611762674840,294.612427581599,294.613630320962,294.614412667777,294.615529087649,294.616354895814,294.617304537353,294.618105317774,294.618999820940,294.619750513456,294.620571877461,294.621303389743,294.622065844177,294.622761441498,294.623494562878,294.624158315837,294.624859266362,294.625511355425,294.626181324809,294.626822290072,294.627474982310,294.628101621750,294.628745269236,294.629360724025,294.629995787723,294.630606676177,294.631233614379,294.631850449854,294.632489276698,294.633131578210,294.633783694722,294.634416923576,294.635060579087,294.635879887091,294.636461305001,294.636884675359,294.637652029189,294.638141139977,294.638841473008,294.639394041814,294.640032897681,294.640570277532,294.641170768377,294.641702798264,294.642261347734,294.642772905792,294.643306345879,294.643794979968,294.644299245108,294.644774186541,294.645252206022,294.645709930009,294.646169809135,294.646611009464,294.647054440013,294.647483844029,294.647911686354,294.648331778385,294.648747546945,294.649158384680,294.649565867274,294.649968353340,294.650369039219,294.651064935732,294.663981094751,294.683925443248,294.721979943793,294.775405736843,294.899963640414,295.052123802178,295.053300501872,294.761550031746,293.769325313208,291.895949233668,289.904069418520,287.578725725334,285.080046474651,280.975685792539];
% INP.AWD = []; %%% when interpolation of input is used

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