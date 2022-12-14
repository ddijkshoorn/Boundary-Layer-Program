%% Experimental data Case_name,'ADA045367_7302_Winter1973_P19'
% Data from ADA045367 case 7302 Profile 19
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
% Low total pressure. Recommende to substitute this case with profile 41
% (Mach 2.2; T0=18.8 Celsius; Re=1.212e6 [1\feet]; cf=149e-5;

%% Data related to extraction of pressure history from drawings
EXP.SP = 451; % measured at station SP, also last station of Mach-number data

%% BL Edge and Wall data (table 1)
EXP.Pn = 19;                % Profile-number
EXP.MaE = 2.1972;           % [-], Mach-number
EXP.ReE = 1.222e6/0.3048;   % [1/m], Unit Reynolds-number
EXP.T0 = 18.5 + 273.15;     % [K]
EXP.cf = 149e-5;            % [-], local skin friction coefficient
EXP.S = 110.4;              % [K], Sutherland's constant for air
EXP.gamma = 1.4;            % [-], specific heat ratio
EXP.R = 287;                % [J/kg/K]

% Obtained from Fernholz & Finley (1977) (they assumed a recovery factor of
% 0.896):
% NB Estimated values!!! Not measured.
EXP.TE = 148.38;            % [K], BL Edge Temperature
EXP.TW = 275.89;            % [K], equals T_recovery (source: Agard)
EXP.P0 = 3.4522e4;          % [Pa], estimated by Fernholz & Finley

%% Profile data (table 1)
% Y-coordinate:
EXP.y = ... % [m]
     0.0254*[0.012,0.014,0.033,0.047,0.057,0.063,0.086,0.100,0.122,0.136,0.168,0.199,0.255,0.295,0.350,0.402,0.443,0.484,0.592,0.689,0.781,0.890,0.992,1.190,1.400,1.597,1.791,1.997,2.192,2.410,2.606,2.786,2.993,3.493,4.492,4.991,5.494];

% Profile data obtained from measured (static) temperature:
EXP.u_ue = ... % [-]
            [0.4870	0.4780	0.5745	0.5957	0.6122	0.6118	0.6395	0.6381	0.657	0.6632	0.6793	0.6915	0.7108	0.7257	0.7376	0.7489	0.7575	0.764	0.7804	0.7911	0.8019	0.8151	0.8241	0.8438	0.8633	0.8783	0.8926	0.9069	0.9203	0.9337	0.9446	0.9552	0.9654	0.9851	0.9984	0.9993	0.9999];
EXP.rho_rhoe =  ... % [-]
            [0.6020	0.5997	0.6050	0.6127	0.6195	0.6193	0.6315	0.658	0.666	0.6705	0.6785	0.685	0.6962	0.7053	0.713	0.72	0.726	0.7302	0.7425	0.7505	0.758	0.769	0.7775	0.7947	0.8132	0.8285	0.845	0.861	0.8772	0.8928	0.9065	0.92	0.9345	0.9653	0.9985	1.0000	1.0000];

% Profile data obtained from measured total temperature:
EXP.Ma = ... % [-]
            [0.8346,0.8133,0.9819,1.0245,1.0587,1.0578,1.1166,1.1373,1.1780,1.1933,1.2294,1.2574,1.3030,1.3391,1.3685,1.3963,1.4181,1.4345,1.4775,1.5059,1.5339,1.5706,1.5966,1.6527,1.7104,1.7565,1.8028,1.8489,1.8939,1.9385,1.9761,2.0130,2.0505,2.1265,2.1921,2.1957,2.1970];

%% Integral properties (table 2c)
EXP.delta_ast = 0.0254*0.9682; % [m], displacement thickness (delta1)
EXP.theta = 0.0254*0.2946; % [m], displacement thickness (delta2)
EXP.delta3 = 0.0254*0.5389; % [m], displacement thickness
EXP.Delta = 0.0254*-0.0128; % [m], enthalpy thickness: integrate this term over BL: (rho u)/(rhoE uE)*(T0/T0E - 1)
EXP.H12 = 3.287; % [-]
EXP.H32 = 1.829; % [-]

% EXP.MaE = 2.1972; % [-], Mach-number, same as table 1
% EXP.ReE = 1.018e5/0.0254; % [1/m], Unit Reynolds-number, same as in table
% 1, but different units!
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%