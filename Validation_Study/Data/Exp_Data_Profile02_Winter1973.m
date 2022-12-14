%% Experimental data Case_name,'ADA045367_7302_Winter1973_P02'
% Data from ADA045367 case 7302 Profile 02
% Original sources (used as much as possible): Winter & Gaudet (1973)
% Composed by D.Dijkshoorn (17-04-2020)

%% Description of case
% Objective was to find how to predict compressible
% boundary layers using incompressible prediction methods corrected with
% compressibility factors. This case is a nice comparison for MM, showing
% the relative compressibility effects inside the BL and in the core (edge)
% flow.
% Winter & Gaudet measured total pressure and total temperature seperately at same single
% station, just before nozzle outlet, at ZPG conditions. cf was also
% measured, T0 and Mach were varied to obtaind data for different unit
% Reynolds-numbers and Mach-numbers.

%% Data obtained from original source
% Mach-number distribution (pressure history) obtained from combining Fig. 6 with Fig. 1 in Winter & Gaudet (1973)
% profile data and integral properties obtained from original source:
% tables 1 (edge data and 'measured temperature' profile data), and 2c.
% The Mach-number obtained from the graph was scaled to make it reach Mach
% 2.8, etc.; The bypass air slots are said to be used above mach 2.2 (Ma > 2.2),
% and thus the simulation of profile 26 (Ma 2.8) might be the most accurate in terms
% of surface length simulated (the other cases might have a longer but unknown surface length, extending upstream).
% For Mach 0.2 they did not list the temperature measurements, since u/ue
% and rho/rhoe were indistuingishable from the tabulated values using
% r=0.89. Thus, this profile data is used here too.
% Furthermore, they say that the later simulations of mach 2.2 are the most
% reliable in their opinion. None of them was used here. I chose to use the
% data with increments of Mach-number for consistency.
% Since pressure and temperature were not measured at the same
% time/experiment, values of cf and P and T were interpolated on basis of
% unit Reynolds-number and location resepctively.

%% Data obtained from AGARD report
%
% Fernholz and Finley consider the temperature measurements as more
% consistent, and thus list the temperature measurements only (pressure
% measurements resulting in Mach-profiles were omitted by them). Measured
% temperatures resulted in listed density and velocity profiles, which were
% used to calculate the other properties. Their aim was to obtain a
% complete data set, I tried to use as much as possible the tabulated data
% from Winter & Gaudet, to avoid having their assumptions affect the
% comparison (validation).

% NB Fernholz & Finley assume a recovery factor. This is not useful for the
% current objective of validation, since it was not at all measured.

%% Doubts
% Note that entry 5 and 6 of Mach and u_ue and rho_rhoe are the same

%% Data related to extraction of pressure history from drawings
EXP.SP = 256; % measured at station SP, also last station of Mach-number data

%% BL Edge and Wall data (table 1)
EXP.Pn = 2;                 % Profile-number
EXP.MaE = 0.2007;           % [-], Mach-number
EXP.ReE = 2.067e6/0.3048;   % [1/m], Unit Reynolds-number
EXP.T0 = 10.8 + 273.15;     % [K]
EXP.cf = 177e-5;            % [-], local skin friction coefficient
EXP.S = 110.4;              % [K], Sutherland's constant for air
EXP.gamma = 1.4;            % [-], specific heat ratio
EXP.R = 287;                % [J/kg/K]

% Obtained from Fernholz & Finley (1977) (they assumed a recovery factor of
% 0.896):
% NB Estimated values!!! Not measured.
EXP.TE = 281.68;            % [K], BL Edge Temperature
EXP.TW = 283.70;            % [K], equals T_recovery (source: Agard)
EXP.P0 = 1.4688e5;          % [Pa], estimated by Fernholz & Finley

%% Profile data (r=0.89) (table 1)
% Y-coordinate:
EXP.y = ... % [m]
     0.0254*[0.012,0.014,0.033,0.047,0.057,0.063,0.086,0.100,0.122,0.136,0.168,0.199,0.255,0.350,0.402,0.443,0.484,0.592,0.689,0.781,0.890,1.190,1.400,1.597,1.791,1.997,2.192,2.410,2.606,2.786,2.993,3.493,3.992,4.492,4.991,5.494,5.991,6.492,6.991];

% Profile data obtained from r=0.89 (not distinguisable from measured (static) temperature according to authors):
EXP.u_ue = ... % [-]
            [0.4914	0.4877	0.5573	0.5734	0.5828	0.5828	0.6190	0.6220	0.6392	0.6532	0.6587	0.6830	0.6935	0.7115	0.7340	0.7315	0.7340	0.7557	0.7581	0.7699	0.7769	0.8086	0.8197	0.8305	0.8455	0.8560	0.8664	0.8767	0.8889	0.8949	0.9048	0.9282	0.9491	0.9639	0.9768	0.9876	0.9966	0.9984	1.0001];
EXP.rho_rhoe =  ... % [-]
            [0.9946	0.9946	0.9951	0.9952	0.9953	0.9953	0.9956	0.9956	0.9958	0.9959	0.9960	0.9962	0.9963	0.9965	0.9967	0.9967	0.9967	0.9969	0.9970	0.9971	0.9972	0.9975	0.9977	0.9978	0.9980	0.9981	0.9982	0.9983	0.9985	0.9986	0.9987	0.9990	0.9993	0.9995	0.9997	0.9998	1.0000	1.0000	1.0000];

% Profile data obtained from measured total temperature:
EXP.Ma = ... % [-]
            [0.0984,0.0976,0.1116,0.1148,0.1167,0.1167,0.1240,0.1246,0.1280,0.1308,0.1319,0.1368,0.1389,0.1425,0.1471,0.1466,0.1471,0.1514,0.1519,0.1543,0.1557,0.1621,0.1643,0.1665,0.1695,0.1716,0.1737,0.1758,0.1783,0.1795,0.1815,0.1862,0.1904,0.1934,0.1960,0.1982,0.2000,0.2004,0.2007];

%% Integral properties (r=0.89) (table 2c)
EXP.delta_ast = 0.0254*0.7071; % [m], displacement thickness (delta1)
EXP.theta = 0.0254*0.5584; % [m], displacement thickness (delta2)
EXP.delta3 = 0.0254*1.0166; % [m], displacement thickness
% EXP.Delta = 0.0254*; % [m], enthalpy thickness: integrate this term over BL: (rho u)/(rhoE uE)*(T0/T0E - 1)
EXP.H12 = 1.266; % [-]
EXP.H32 = 1.820; % [-]

% EXP.MaE = 0.2007; % [-], Mach-number, same as table 1
% EXP.ReE = 1.722e5/0.0254; % [1/m], Unit Reynolds-number, same as in table
% 1, but different units!
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%