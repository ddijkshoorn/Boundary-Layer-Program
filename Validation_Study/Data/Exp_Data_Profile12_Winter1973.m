%% Experimental data Case_name,'ADA045367_7302_Winter1973_P12'
% Data from ADA045367 case 7302 Profile 12
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
% 

%% Data related to extraction of pressure history from drawings
EXP.SP = 396; % measured at station SP, also last station of Mach-number data

%% BL Edge and Wall data (table 1)
EXP.Pn = 12;                % Profile-number
EXP.MaE = 1.4003;           % [-], Mach-number
EXP.ReE = 2.043e6/0.3048;   % [1/m], Unit Reynolds-number
EXP.T0 = 18.6 + 273.15;     % [K]
EXP.cf = 165e-5;            % [-], local skin friction coefficient
EXP.S = 110.4;              % [K], Sutherland's constant for air
EXP.gamma = 1.4;            % [-], specific heat ratio
EXP.R = 287;                % [J/kg/K]

% Obtained from Fernholz & Finley (1977) (they assumed a recovery factor of
% 0.896):
% NB Estimated values!!! Not measured.
EXP.TE = 209.57;            % [K], BL Edge Temperature
EXP.TW = 282.72;            % [K], equals T_recovery (source: Agard)
EXP.P0 = 9.5936e4;          % [Pa], estimated by Fernholz & Finley

%% Profile data (table 1)
% Y-coordinate:
EXP.y = ... % [m]
     0.0254*[0.012,0.014,0.033,0.047,0.057,0.063,0.086,0.100,0.122,0.136,0.168,0.199,0.255,0.295,0.350,0.400,0.443,0.484,0.592,0.689,0.781,0.890,0.992,1.190,1.400,1.597,1.791,1.997,2.192,2.410,2.606,2.786,2.993,3.493,3.992,4.492,4.991,5.494];

% Profile data obtained from measured (static) temperature:
EXP.u_ue = ... % [-]
            [0.5242	0.5172	0.5654	0.5868	0.5974	0.6021	0.6288	0.6395	0.6555	0.6682	0.6815	0.6964	0.7157	0.7271	0.7425	0.7511	0.7584	0.7647	0.7832	0.7941	0.8039	0.8132	0.8261	0.8439	0.8621	0.8764	0.8906	0.9040	0.9167	0.9289	0.9388	0.9486	0.9583	0.9784	0.9941	0.9993	1.0004	1.0000];

EXP.rho_rhoe =  ... % [-]
            [0.7950	0.7945	0.8053	0.8104	0.8145	0.8150	0.8225	0.8259	0.8305	0.8350	0.8397	0.8435	0.8504	0.8550	0.8601	0.8635	0.8660	0.8695	0.8765	0.8810	0.8853	0.8890	0.8953	0.9040	0.9125	0.9199	0.9267	0.9345	0.9408	0.9482	0.9539	0.9598	0.9653	0.9782	0.9875	0.9939	0.9974	1.0000];

% Profile data obtained from measured total temperature:
EXP.Ma = ... % [-]
            [0.6544,0.6456,0.7104,0.7397,0.7549,0.7611,0.7985,0.8139,0.8365,0.8550,0.8745,0.8955,0.9242,0.9415,0.9643,0.9774,0.9882,0.9984,1.0268,1.0437,1.0592,1.0737,1.0946,1.1235,1.1531,1.1770,1.2005,1.2237,1.2450,1.2666,1.2840,1.3014,1.3184,1.3551,1.3832,1.3951,1.3990,1.4003];

%% Integral properties (table 2c)
EXP.delta_ast = 0.0254*0.7278; % [m], displacement thickness (delta1)
EXP.theta = 0.0254*0.3431; % [m], displacement thickness (delta2)
EXP.delta3 = 0.0254*0.6264; % [m], displacement thickness
EXP.Delta = 0.0254*0.0010; % [m], enthalpy thickness: integrate this term over BL: (rho u)/(rhoE uE)*(T0/T0E - 1)
EXP.H12 = 2.121; % [-]
EXP.H32 = 1.826; % [-]

% EXP.MaE = 1.4003; % [-], Mach-number, same as table 1
% EXP.ReE = 1.702e5/0.0254; % [1/m], Unit Reynolds-number, same as in table
% 1, but different units!
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%